#include <votca/ctp/ewaldnd.h>
#include <votca/ctp/poissongrid.h>
#include <boost/format.hpp>
#include <algorithm>
#include <boost/math/special_functions/round.hpp>
#include <boost/timer/timer.hpp>


namespace votca { namespace ctp {

using boost::format;

    

Ewald3DnD::~Ewald3DnD() {
    vector< PolarSeg* >::iterator sit;
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit)
        delete (*sit);
    
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();
    
    _polar_qm0.clear();
    _polar_mm1.clear();
    _polar_mm2.clear();
    
    delete _fg_table;
    _fg_table = 0;
}
    
    
Ewald3DnD::Ewald3DnD(Topology *top, PolarTop *ptop, Property *opt, Logger *log) 
    : _top(top), _ptop(ptop), _log(log), _fg_table(0) {
    
    // EVALUATE OPTIONS
    string pfx = "options.ewald";
    // Ewald parameters
    string cmethod = opt->get(pfx+".coulombmethod.method").as<string>();
    assert(cmethod == "ewald" && "<::Ewald3DnD> CMETHOD NOT IMPLEMENTED");
    if (opt->exists(pfx+".coulombmethod.cutoff")) {
        _R_co = opt->get(pfx+".coulombmethod.cutoff").as<double>();
    }
    else
        _R_co = 6.;
    if (opt->exists(pfx+".coulombmethod.shape")) {
        _shape = opt->get(pfx+".coulombmethod.shape").as<string>();
    }
    else
        _shape = "xyslab";
    if (opt->exists(pfx+".coulombmethod.save_nblist")) {
        _save_nblist = opt->get(pfx+".coulombmethod.save_nblist").as<bool>();
    }
    else
        _save_nblist = true;
    // Convergence
    if (opt->exists(pfx+".convergence.energy")) {
        _crit_dE = opt->get(pfx+".convergence.energy").as<double>();
    }
    else
        _crit_dE = 1e-5;
    if (opt->exists(pfx+".convergence.kfactor"))
        _kfactor = opt->get(pfx+".convergence.kfactor").as<double>();
    else
        _kfactor = 100.;
    if (opt->exists(pfx+".convergence.rfactor"))
        _rfactor = opt->get(pfx+".convergence.rfactor").as<double>();
    else
        _rfactor = 6.;
    // Polar parameters
    string pmethod = opt->get(pfx+".coulombmethod.method").as<string>();
    assert(pmethod == "ewald" && "<::Ewald3DnD> PMETHOD NOT IMPLEMENTED");
    if (opt->exists(pfx+".polarmethod.induce"))
        _polar_do_induce = opt->get(pfx+".polarmethod.induce").as<bool>();
    else
        _polar_do_induce = false;
    if (opt->exists(pfx+".polarmethod.cutoff")) 
        _polar_cutoff = opt->get(pfx+".polarmethod.cutoff").as<double>();
    else
        _polar_cutoff = 0.0;
    if (opt->exists(pfx+".polarmethod.wSOR_N"))
        _polar_wSOR_N = opt->get(pfx+".polarmethod.wSOR_N").as<double>();
    else
        _polar_wSOR_N = 0.35;
    if (opt->exists(pfx+".polarmethod.wSOR_C"))
        _polar_wSOR_C = opt->get(pfx+".polarmethod.wSOR_C").as<double>();
    else
        _polar_wSOR_C = 0.30;
    if (opt->exists(pfx+".polarmethod.aDamp"))
        _polar_aDamp = opt->get(pfx+".polarmethod.aDamp").as<double>();
    else
        _polar_aDamp = 0.390;
    // Coarse-graining
    if (opt->exists(pfx+".coarsegrain.cg_background")) {
        _coarse_do_cg_background = 
            opt->get(pfx+".coarsegrain.cg_background").as<bool>();
    }
    else
        _coarse_do_cg_background = false;
    if (opt->exists(pfx+".coarsegrain.cg_foreground")) {
        _coarse_do_cg_foreground =
            opt->get(pfx+".coarsegrain.cg_foreground").as<bool>();
    }
    else
        _coarse_do_cg_foreground = false;
    if (opt->exists(pfx+".coarsegrain.cg_radius")) {
        _coarse_cg_radius =
            opt->get(pfx+".coarsegrain.cg_radius").as<double>();
    }
    else
        _coarse_cg_radius = _polar_cutoff;
    if (opt->exists(pfx+".coarsegrain.cg_anisotropic")) {
        _coarse_cg_anisotropic =
            opt->get(pfx+".coarsegrain.cg_anisotropic").as<bool>();
    }
    else
        _coarse_cg_anisotropic = false;
    // Tasks to perform
    if (opt->exists(pfx+".tasks.calculate_fields")) {
        _task_calculate_fields 
            = opt->get(pfx+".tasks.calculate_fields").as<bool>();
    }
    else
        _task_calculate_fields = false;
    if (opt->exists(pfx+".tasks.polarize_fg")) {
        _task_polarize_fg = opt->get(pfx+".tasks.polarize_fg").as<bool>();
    }
    else
        _task_polarize_fg = false;
    if (opt->exists(pfx+".tasks.evaluate_energy")) {
        _task_evaluate_energy
            = opt->get(pfx+".tasks.evaluate_energy").as<bool>();
    }
    else
        _task_evaluate_energy = false;
    if (opt->exists(pfx+".tasks.solve_poisson")) {
        _task_solve_poisson
            = opt->get(pfx+".tasks.solve_poisson").as<bool>();
    }
    else
        _task_solve_poisson = false;
    
    // EWALD INTERACTION PARAMETERS (GUESS ONLY)
    _K_co = _kfactor/_R_co;
    _alpha = _rfactor/_R_co;
    _ewdactor = EwdInteractor(_alpha, _polar_aDamp);
    
    _did_field_pin_R_shell = false;
    _did_generate_kvectors = false;
    
    // SET-UP REAL & RECIPROCAL SPACE
    _a = _top->getBox().getCol(0);
    _b = _top->getBox().getCol(1);
    _c = _top->getBox().getCol(2);
    _LxLyLz = _a*(_b^_c);
    _LxLy = abs(_a ^ _b);
    
    _A = 2*M_PI/_LxLyLz * _b^_c;
    _B = 2*M_PI/_LxLyLz * _c^_a;
    _C = 2*M_PI/_LxLyLz * _a^_b;

    _na_max = ceil((_R_co+_polar_cutoff)/maxnorm(_a)-0.5)+1;
    _nb_max = ceil((_R_co+_polar_cutoff)/maxnorm(_b)-0.5)+1;
    _nc_max = ceil((_R_co+_polar_cutoff)/maxnorm(_c)-0.5)+1;

    _NA_max = ceil(_K_co/maxnorm(_A));
    _NB_max = ceil(_K_co/maxnorm(_B));
    _NC_max = ceil(_K_co/maxnorm(_C));
    
    // SET-UP POLAR GROUNDS (FORE-, MID-, BACK-)
    _center = ptop->getCenter();
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();
    
    _fg_C = ptop->FGC();
    _fg_N = ptop->FGN();
    _bg_N = ptop->BGN();
    
    // Construct periodic neutral background
    _bg_P.insert(_bg_P.end(), _fg_N.begin(), _fg_N.end());
    _bg_P.insert(_bg_P.end(), _bg_N.begin(), _bg_N.end());    
    // Grow foreground according to induction cut-off
    this->ExpandForegroundReduceBackground(_polar_cutoff);
    // Coarse-grain as demanded by input
    boost::timer::cpu_timer cpu_t;
    cpu_t.start();
    boost::timer::cpu_times t0 = cpu_t.elapsed();
    this->CoarseGrainDensities(_coarse_do_cg_background, 
        _coarse_do_cg_foreground, _coarse_cg_radius);
    boost::timer::cpu_times t1 = cpu_t.elapsed();
    _t_coarsegrain = (t1.wall-t0.wall)/1e9/60.;
    
    // SET-UP MIDGROUND (INCLUDING PERIODIC IMAGES IF REQUIRED)
    LOG(logINFO,*_log) << flush;
    LOG(logINFO,*_log) << "Generate periodic images. ";
    this->SetupMidground(_R_co);
    
    // CALCULATE COG POSITIONS, NET CHARGE
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    double Q_fg_C = 0.0;
    double Q_fg_C_2nd = 0.0;
    double Q_fg_N = 0.0;
    double Q_mg_N = 0.0;
    double Q_bg_N = 0.0;
    double Q_bg_P = 0.0;  
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        (*sit)->CalcPos();
        double Qseg = (*sit)->CalcTotQ();
        Q_fg_C += Qseg;
        Q_fg_C_2nd += Qseg*Qseg / _fg_C.size();
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_fg_N += (*sit)->CalcTotQ();
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_mg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_P += (*sit)->CalcTotQ();
    }
    
    // DETERMINE JOB TYPE
    int iQ1 = boost::math::iround(Q_fg_C);
    int iQ2 = boost::math::iround(Q_fg_C_2nd);
    if (iQ1 == +1)         _jobType = "hole-like";
    else if (iQ1 == -1)    _jobType = "electron-like";
    else if (iQ1 == 0 && iQ2 == 0) _jobType = "neutral";
    else if (iQ1 == 0 && iQ2 > 0) _jobType = "charge-transfer-like";
    else _jobType = "bipolaron-like";

    
    LOG(logINFO,*_log)
        << (format("Net ground charge and size:")).str()
        << flush << (format("  o Q(FGC) = %1$+1.3fe |FGC| = %2$+5d") % Q_fg_C % _fg_C.size()).str()
        << flush << (format("  o Q(FGN) = %1$+1.3fe |FGN| = %2$+5d") % Q_fg_N % _fg_N.size()).str()
        << flush << (format("  o Q(MGN) = %1$+1.3fe |MGN| ~ %2$+5d") % Q_mg_N % _mg_N.size()).str()
        << flush << (format("  o Q(BGN) = %1$+1.3fe |BGN| = %2$+5d") % Q_bg_N % _bg_N.size()).str()
        << flush << (format("  o Q(BGP) = %1$+1.3fe |BGP| = %2$+5d") % Q_bg_P % _bg_P.size()).str()
        << flush << (format("  o Job type '%3$s' (iQ1=%1$d, iQ2=%2$d)") % iQ1 % iQ2 % _jobType).str()
        << flush;
    
    if (std::abs(Q_bg_P) > 1e-2) {
        cout << endl;
        cout << endl << format("***************************** ERROR ******************************");
        cout << endl << format("       Background charge |Q(BGP)| is larger than 0.01e.");
        cout << endl << format("       Be more precise: e.g. rounding error?");
        cout << endl << format("       Or think again: e.g. erroneous parametrization?");
        cout << endl << format("******************************************************************");
        cout << endl;
    }
    
    // CALCULATE NET DIPOLE OF BGP & FGC
    vec netdpl_bgP = vec(0,0,0);
    vec netdpl_fgC = vec(0,0,0);
    double qzz_bgP = 0.0;
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl_bgP += (*pit)->getPos() * (*pit)->getQ00();
            qzz_bgP += (*pit)->getQ00() * ((*pit)->getPos().getZ() * (*pit)->getPos().getZ());
        }
    }
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl_fgC += (*pit)->getPos() * (*pit)->getQ00();
        }
    }
    
    LOG(logINFO,*_log)
        << (format("Net dipole moment of background density")).str()
        << flush << (format("  o D(BGP) [e*nm]           = %1$+1.3f %2$+1.3f %3$+1.3f  ") 
        % netdpl_bgP.getX() % netdpl_bgP.getY() % netdpl_bgP.getZ()).str();
    LOG(logINFO,*_log)
        << flush << (format("  o D(FGC) [e*nm]           = %1$+1.3f %2$+1.3f %3$+1.3f  ") 
        % netdpl_fgC.getX() % netdpl_fgC.getY() % netdpl_fgC.getZ()).str();
    LOG(logINFO,*_log)
        << flush << (format("  o Sigma q|z|**2 [e*nm**2] = %1$+1.7f   ")
        % qzz_bgP) << flush;
    
    
    // ZERO ENERGIES
    _ER  = EWD::triple<>(0, 0, 0);
    _EC  = EWD::triple<>(0, 0, 0);
    _EK  = EWD::triple<>(0, 0, 0);
    _E0  = EWD::triple<>(0, 0, 0);
    _ET  = EWD::triple<>(0, 0, 0);
    _EDQ = EWD::triple<>(0, 0, 0);
    _EJ  = EWD::triple<>(0, 0, 0);
    _polar_ETT = 0;
    _polar_EPP = 0;  _polar_EPU = 0;  _polar_EUU = 0;
    _polar_EF00 = 0; _polar_EF01 = 0; _polar_EF02 = 0;
    _polar_EF11 = 0; _polar_EF12 = 0;
    _polar_EM0 = 0;  _polar_EM1 = 0;  _polar_EM2 = 0;
    _Estat = 0;      _Eindu = 0;      _Eppuu = 0;
    
    return;
}


void Ewald3DnD::ExpandForegroundReduceBackground(double polar_R_co) {
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "Set-up polar grounds" << flush;
    
    vector<PolarSeg*>::iterator sit1;
    vector<PolarSeg*>::iterator sit2;
    
    assert(_polar_mm1.size() == 0 
        && _polar_qm0.size() == 0 
        && _polar_mm2.size() == 0);
    
    // Target containers
    vector<PolarSeg*> exp_fg_C;
    vector<PolarSeg*> exp_fg_N;
    vector<PolarSeg*> red_bg_N;
    
    // Image boxes to consider, set-up boolean foreground table
    int polar_na_max = ceil(polar_R_co/maxnorm(_a)-0.5)+1;
    int polar_nb_max = ceil(polar_R_co/maxnorm(_b)-0.5)+1;
    int polar_nc_max = ceil(polar_R_co/maxnorm(_c)-0.5)+1;
    
    LOG(logDEBUG,*_log) << "  o Expanding cell space for neighbour search:"
            " +/-" << polar_na_max << " x +/-" << polar_nb_max << " x +/-" 
            << polar_nc_max << flush;
    
    _fg_table = new ForegroundTable(
        _bg_P.size(), polar_na_max, polar_nb_max, polar_nc_max);
    
    // Max. distance between any two segments in FGC before expansion
    _max_int_dist_qm0 = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = sit1+1; sit2 < _fg_C.end(); ++sit2) {
            double R = votca::tools::abs((*sit1)->getPos()-(*sit2)->getPos());
            _max_int_dist_qm0 = (R > _max_int_dist_qm0) ? R : _max_int_dist_qm0;
        }
    }
    
    LOG(logDEBUG,*_log) << "  o Max. int. distance (QM0): " << 
        _max_int_dist_qm0 << "nm" << flush;
    
    // Foreground remains in foreground + contributes to QM0
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        exp_fg_C.push_back(*sit1);
        _polar_qm0.push_back(*sit1);
        _fg_table->AddToForeground((*sit1)->getId(), 0, 0, 0);
    }
    for (sit1 = _fg_N.begin(); sit1 < _fg_N.end(); ++sit1) {
        exp_fg_N.push_back(*sit1);
    }    
    
    // Background expands and migrates to foreground OR remains background
    int allocated_count_n = 0;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        PolarSeg *seg_bg = *sit1;
        for (int na = -polar_na_max; na < polar_na_max+1; ++na) {
        for (int nb = -polar_nb_max; nb < polar_nb_max+1; ++nb) {
        for (int nc = -polar_nc_max; nc < polar_nc_max+1; ++nc) {            
            vec L = na*_a + nb*_b + nc*_c;
            
            bool in_central_cell = (na == 0 && nb == 0 && nc == 0);
            bool identical = false;
            bool within_range = false;
            
            for (sit2 = _fg_C.begin(); sit2 < _fg_C.end(); ++sit2) {
                PolarSeg *seg_fg = *sit2;
                // Identical ?
                if (in_central_cell && seg_bg->getId() == seg_fg->getId()) {
                    assert(identical==false);
                    identical = true;
                }
                // Within range ?
                // NOTE We have to truncate the decimal places for the radius
                // and the cut-off before we draw the comparison - otherwise
                // numerics may play a prank on us and produce different
                // foregrounds for different charge states, making energy
                // differences for those states pretty much meaningless.
                // Note that abs(dR_L-polar_R_co) < 1e-xy won't do here,
                // since we then to a large degree still rely on machine 
                // precision.
                double dR_L = votca::tools::abs(
                    seg_bg->getPos() + L - seg_fg->getPos());
                // Compare up to 3 decimal places (= 1e-3 nm)
                if (int(dR_L*1e3+0.5) <= int(polar_R_co*1e3+0.5)) {
                    within_range = true;
                }
            }
            
            if (!identical && within_range) {
                _fg_table->AddToForeground(seg_bg->getId(), na, nb, nc);
                // Add new shifted clone to fgC and MM1, depolarize = true
                PolarSeg *seg_bg_clone_fgc = new PolarSeg(seg_bg, true);
                seg_bg_clone_fgc->Translate(L);
                exp_fg_C.push_back(seg_bg_clone_fgc);
                _polar_mm1.push_back(seg_bg_clone_fgc);
                // Add original to fgN OR a shifted clone if image box != 0,
                // depolarize = false
                if (in_central_cell)
                    exp_fg_N.push_back(seg_bg);
                else {
                    allocated_count_n += 1;
                    PolarSeg *seg_bg_clone_fgn = new PolarSeg(seg_bg, false);
                    seg_bg_clone_fgn->Translate(L);
                    exp_fg_N.push_back(seg_bg_clone_fgn);
                }
            }
            else if (!identical && !within_range) {
                // Add original to bgN OR a shifted clone if image box != 0
                if (in_central_cell)
                    red_bg_N.push_back(seg_bg);
                else {
                    ;
                    //PolarSeg *seg_bg_clone_bgn = new PolarSeg(seg_bg);
                    //seg_bg_clone_bgn->Translate(L);
                    //red_bg_N.push_back(seg_bg_clone_bgn);
                }
            }
            else {
                ;
            }
        }}} // Sum over image boxes
    } // Sum over periodic neutral background
    
    // Exchange new for old containers
    _fg_C.clear();
    _fg_N.clear();
    _bg_N.clear();
    _fg_C = exp_fg_C;
    _fg_N = exp_fg_N;
    _bg_N = red_bg_N;
    
    bool clean = true;
    _ptop->setFGC(_fg_C, true);
    _ptop->setFGN(_fg_N, true);
    _ptop->setBGN(_bg_N, true);
    clean = false; // Already deleted via fgC
    _ptop->setQM0(_polar_qm0, clean);
    _ptop->setMM1(_polar_mm1, clean);
    _ptop->setMM2(_polar_mm2, clean);    
    
    // Sanity checks
    assert(_polar_qm0.size()+_polar_mm1.size() == _fg_C.size());
    assert(_fg_C.size() == _fg_N.size());
    
    return;
}


void Ewald3DnD::CoarseGrainDensities(bool cg_bg, bool cg_fg, double cg_radius) {
    
    LOG(logDEBUG,*_log) << "Coarse-graining agenda" << flush;
    LOG(logDEBUG,*_log) << "  o Coarse-grain background:   " << ((_coarse_do_cg_background) ? "yes" : "no") << flush;
    LOG(logDEBUG,*_log) << "  o Coarse-grain foreground:   " << ((_coarse_do_cg_foreground) ? "yes" : "no") << flush;
    LOG(logDEBUG,*_log) << "  o Anisotropic P-tensor:      " << ((_coarse_cg_anisotropic) ? "yes" : "no") << flush;
    LOG(logDEBUG,*_log) << "  o Coarse-graining radius:    " << _coarse_cg_radius << " nm" << flush;
    
    // COARSE-GRAIN BACKGROUND
    if (cg_bg) {
        LOG(logDEBUG,*_log) << "Coarse-grain background" << flush;
        
        int count_bgp = 0;
        int count_fgn = 0;
        int count_fgc = 0;
        int count_bgp_id_in_fg = 0;

        for (vector<PolarSeg*>::iterator sit = _bg_P.begin();
            sit != _bg_P.end(); ++sit) {
            if (_fg_table->IsInForeground((*sit)->getId(), 0, 0, 0)) {
                ++count_bgp_id_in_fg;
            }
            // By checking whether there are more polar sites than polar
            // fragments in the segment, one avoids coarse-graining twice
            assert((*sit)->size() >= (*sit)->PolarFrags().size()
               && "<::CoarseGrainDensities> BGN: FEWER POLAR SITES THAN FRAGS");
            if ((*sit)->size() > (*sit)->PolarFrags().size()) {
                //cout << "\rMST DBG ...   o BGP ID = " << (*sit)->getId() 
                //     << "   " << flush;
                (*sit)->Coarsegrain(_coarse_cg_anisotropic);
                ++count_bgp;
            }
        }

        for (vector<PolarSeg*>::iterator sit = _fg_N.begin();
            sit != _fg_N.end(); ++sit) {
            // By checking whether there are more polar sites than polar
            // fragments in the segment, one avoids coarse-graining twice
            assert((*sit)->size() >= (*sit)->PolarFrags().size()
               && "<::CoarseGrainDensities> FGN: FEWER POLAR SITES THAN FRAGS");
            if ((*sit)->size() > (*sit)->PolarFrags().size()) {
                //cout << "\rMST DBG ...   o FGN ID = " << (*sit)->getId() 
                //     << "   " << flush;
                (*sit)->Coarsegrain(_coarse_cg_anisotropic);
                ++count_fgn;
            }
        }
        
        LOG(logDEBUG,*_log) << "  o Coarse-grained "
             << count_bgp << "/" << _bg_P.size() << " (BGP) "
             << count_fgn << "/" << _fg_N.size() << " (FGN) "
             << "with " << count_bgp_id_in_fg << " FG-table counts" << flush;
    }
    // COARSE-GRAIN FOREGROUND
    if (cg_fg) {
        LOG(logDEBUG,*_log) << "Coarse-grain foreground" << flush;
        
        int count_fgc = 0;
        
        for (vector<PolarSeg*>::iterator sit = _fg_C.begin();
            sit != _fg_C.end(); ++sit) {
            assert((*sit)->size() >= (*sit)->PolarFrags().size()
               && "<::CoarseGrainDensities> FGC: FEWER POLAR SITES THAN FRAGS");            
            // Only coarse-grain the site if the minimum distance of approach
            // to any of the sites in QM0 is larger than cg_radius
            // Negative radius: Coarse-grain all sites in FGC
            double min_R = std::abs(2*cg_radius);
            for (vector<PolarSeg*>::iterator qit = _polar_qm0.begin();
                qit != _polar_qm0.end(); ++qit) {
                double R = votca::tools::abs((*sit)->getPos()-(*qit)->getPos());
                if (R < min_R) min_R = R;
            }

			if (!(min_R <= _polar_cutoff+1e-3)) {
				cout << endl << "ASSERTION IMMINENT" << endl;
			}
            assert(min_R <= _polar_cutoff+1e-3
                && "<::CoarseGrainDensities> FGC: INCONSISTENT WITH EXPANSION");
            // Different from above, here there should be no danger of 
            // coarse-graining twice
            if (min_R > cg_radius) {
                //cout << "\rMST DBG ...   o FGC ID = " << (*sit)->getId() 
                //     << "   " << flush;
                (*sit)->Coarsegrain(_coarse_cg_anisotropic);
                ++count_fgc;
            }
        }
        
        LOG(logDEBUG,*_log) << "  o Coarse-grained "
             << count_fgc << "/" << _fg_C.size() << " (FGC) " << flush;
    }    
    return;
}


void Ewald3DnD::SetupMidground(double R_co) {
    // SET-UP MIDGROUND
    // NOTE No periodic-boundary correction here: We require that all
    //      polar-segment CoM coords. be folded with respect to the central
    //      segment
    // NOTE Excludes interaction of fg polar segments with neutral selfs in 
    //      real-space sum
    // NOTE Includes periodic images if within cut-off    

    // CLEAR ANY EXTANT MIDGROUND
    vector<PolarSeg*>::iterator sit;
    vector<PolarSeg*>::iterator sit2;
    vector<APolarSite*> ::iterator pit;
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit)
        delete (*sit);
    _mg_N.clear();
    
    // SAMPLE MIDGROUND FROM BGP EXCLUDING CENTRAL SEG.
    assert(_fg_N.size() == _fg_C.size());
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        PolarSeg *pseg = *sit;
        // Periodic images
        for (int na = -_na_max; na < _na_max+1; ++na) {
        for (int nb = -_nb_max; nb < _nb_max+1; ++nb) {
        for (int nc = -_nc_max; nc < _nc_max+1; ++nc) {
            vec L = na*_a + nb*_b + nc*_c;
            // In foreground ?
            bool is_in_fg = _fg_table->IsInForeground(pseg->getId(),na,nb,nc);            
            if (!is_in_fg) {
                bool is_within_range = false;
                // Within range ?
                for (sit2 = _fg_C.begin(); sit2 < _fg_C.end(); ++sit2) {
                    vec pos_L = pseg->getPos() + L;
                    double dR_L = abs((*sit2)->getPos()-pos_L);
                    if (dR_L <= R_co) {
                        is_within_range = true;
                    }
                }
                // Add if appropriate, depolarize = false
                if (is_within_range) {
                    PolarSeg *newSeg = new PolarSeg(pseg, false);
                    newSeg->Translate(L);
                    _mg_N.push_back(newSeg);
                }
            }
            else ;
        }}} // Loop over na, nb, nc
    } // Loop over BGP
    
    return;
}


void Ewald3DnD::WriteDensitiesPDB(string pdbfile) {
    // COORDINATE OUTPUT FOR VISUAL CHECK
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;    
    FILE *out;
    out = fopen(pdbfile.c_str(),"w");
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGC");
        }
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGN");
        }
    }
    for (sit = _polar_qm0.begin(); sit < _polar_qm0.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "QM0");
        }
    }
    for (sit = _polar_mm1.begin(); sit < _polar_mm1.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MM1");
        }
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MGN");
        }
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGN");
        }
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGP");
        }
    }
    fclose(out);
    return;
}


void Ewald3DnD::WriteDensitiesPtop(string fg, string mg, string bg) {
    // FGC, FGN, BGN, QM0, MM1, MM2
    _ptop->SaveToDrive(fg);
    // MGN
    PolarTop mg_ptop;
    mg_ptop.setBGN(_mg_N, false);
    mg_ptop.SaveToDrive(mg);
    // BGP
    PolarTop bg_ptop;
    bg_ptop.setBGN(_bg_P, false);
    bg_ptop.SaveToDrive(bg);
    
    string fg_pdb = fg + ".pdb";
    string mg_pdb = mg + ".pdb";
    string bg_pdb = bg + ".pdb";
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;    
    FILE *out;
    out = fopen(fg_pdb.c_str(),"w");
    for (sit = _polar_qm0.begin(); sit < _polar_qm0.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "QM0");
        }
    }
    for (sit = _polar_mm1.begin(); sit < _polar_mm1.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MM1");
        }
    }
    fclose(out);
    out = fopen(mg_pdb.c_str(),"w");
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MGN");
        }
    }
    fclose(out);
    out = fopen(bg_pdb.c_str(),"w");
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGP");
        }
    }
    fclose(out);
    return;
}

void Ewald3DnD::Evaluate() {
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "System & Ewald parameters (" << IdentifyMethod() << ")" << flush;
    LOG(logDEBUG,*_log) << "  o Real-space unit cell:      " << _a << " x " << _b << " x " << _c << flush;
    LOG(logDEBUG,*_log) << "  o Real-space c/o (guess):    " << _R_co << " nm" << flush;
    LOG(logDEBUG,*_log) << "  o na(max), nb(max), nc(max): " << _na_max << ", " << _nb_max << ", " << _nc_max << flush;
    LOG(logDEBUG,*_log) << "  o 1st Brillouin zone:        " << _A << " x " << _B << " x " << _C << flush;
    LOG(logDEBUG,*_log) << "  o Reciprocal-space c/o:      " << _K_co << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o R-K switching param.       " << _alpha << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o Unit-cell volume:          " << _LxLyLz << " nm**3" << flush;
    LOG(logDEBUG,*_log) << "  o LxLy (for 3D2D EW):        " << _LxLy << " nm**2" << flush;
    LOG(logDEBUG,*_log) << "  o kx(max), ky(max), kz(max): " << _NA_max << ", " << _NB_max << ", " << _NC_max << flush;
    
    LOG(logDEBUG,*_log) << "Tasks to perform (" << IdentifyMethod() << ")" << flush;
    LOG(logDEBUG,*_log) << "  o Calculate fg fields:       " << ((_task_calculate_fields) ? "yes" : "no") << flush;
    LOG(logDEBUG,*_log) << "  o Polarize foreground:       " << ((_task_polarize_fg) ? "yes" : "no") << flush;
    LOG(logDEBUG,*_log) << "  o Evaluate energy:           " << ((_task_evaluate_energy) ? "yes" : "no") << flush;
        
    // TEASER OUTPUT PERMANENT FIELDS
    LOG(logDEBUG,*_log) << flush << "Background fields (BGP):" << flush;
    int fieldCount = 0;
    for (vector<PolarSeg*>::iterator sit1 = _bg_P.begin()+288; sit1 < _bg_P.end(); ++sit1) {
        PolarSeg *pseg = *sit1;
        Segment *seg = _top->getSegment(pseg->getId());
        LOG(logDEBUG,*_log) << "ID = " << pseg->getId() << " (" << seg->getName() << ") " << flush;
        for (PolarSeg::iterator pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
            vec fp = (*pit1)->getFieldP();
            vec fu = (*pit1)->getFieldU();
            vec u1 = (*pit1)->getU1();
            LOG(logDEBUG,*_log)
               << (format("FPU = (%1$+1.7e %2$+1.7e %3$+1.7e) V/m    ") 
                    % (fp.getX()*EWD::int2V_m+fu.getX()*EWD::int2V_m)
                    % (fp.getY()*EWD::int2V_m+fu.getY()*EWD::int2V_m) 
                    % (fp.getZ()*EWD::int2V_m+fu.getZ()*EWD::int2V_m)).str();
            LOG(logDEBUG,*_log)
               << (format("U1* = (%1$+1.7e %2$+1.7e %3$+1.7e) e*nm") 
                    % (u1.getX())
                    % (u1.getY()) 
                    % (u1.getY())).str() << flush;
            fieldCount += 1;
            if (fieldCount > 10) {
                LOG(logDEBUG,*_log)
                    << "FPU = ... ... ..." << flush;
                break;
            }
        }
        if (fieldCount > 10) break;
    }
        
    boost::timer::cpu_timer cpu_t;
    cpu_t.start();
    boost::timer::cpu_times t0 = cpu_t.elapsed();
    if (_task_calculate_fields) EvaluateFields();
    boost::timer::cpu_times t1 = cpu_t.elapsed();
    if (_task_polarize_fg) EvaluateInduction();
    boost::timer::cpu_times t2 = cpu_t.elapsed();
    if (_task_evaluate_energy) EvaluateEnergy();
    boost::timer::cpu_times t3 = cpu_t.elapsed();
    if (_task_solve_poisson) EvaluatePoisson();
    
    _t_fields    = (t1.wall-t0.wall)/1e9/60.;
    _t_induction = (t2.wall-t1.wall)/1e9/60.;
    _t_energy    = (t3.wall-t2.wall)/1e9/60.;
    
    // TEASER OUTPUT PERMANENT FIELDS
    LOG(logDEBUG,*_log) << flush << "Background fields (BGP):" << flush;
    fieldCount = 0;
    for (vector<PolarSeg*>::iterator sit1 = _bg_P.begin()+288; sit1 < _bg_P.end(); ++sit1) {
        PolarSeg *pseg = *sit1;
        Segment *seg = _top->getSegment(pseg->getId());
        LOG(logDEBUG,*_log) << "ID = " << pseg->getId() << " (" << seg->getName() << ") " << flush;
        for (PolarSeg::iterator pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
            vec fp = (*pit1)->getFieldP();
            vec fu = (*pit1)->getFieldU();
            vec u1 = (*pit1)->getU1();
            LOG(logDEBUG,*_log)
               << (format("FPU = (%1$+1.7e %2$+1.7e %3$+1.7e) V/m    ") 
                    % (fp.getX()*EWD::int2V_m+fu.getX()*EWD::int2V_m)
                    % (fp.getY()*EWD::int2V_m+fu.getY()*EWD::int2V_m) 
                    % (fp.getZ()*EWD::int2V_m+fu.getZ()*EWD::int2V_m)).str();
            LOG(logDEBUG,*_log)
               << (format("U1* = (%1$+1.7e %2$+1.7e %3$+1.7e) e*nm") 
                    % (u1.getX())
                    % (u1.getY()) 
                    % (u1.getY())).str() << flush;
            fieldCount += 1;
            if (fieldCount > 10) {
                LOG(logDEBUG,*_log)
                    << "FPU = ... ... ..." << flush;
                break;
            }
        }
        if (fieldCount > 10) break;
    }
    
    
    
    
    double outer_epp = _ET._pp;
    double outer_eppu = _ET._pu + _ET._uu;
    double outer = outer_epp + outer_eppu;
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logINFO,*_log)
        << (format("Interaction FGC -> ***")).str()
        << flush << (format("  + EPP(FGC->MGN)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _ER.Sum() % _ER._pp % _ER._pu % _ER._uu).str()
        << flush << (format("  + EKK(FGC->BGP)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _EK.Sum() % _EK._pp % _EK._pu % _EK._uu).str()       
        << flush << (format("  - EPP(FGC->FGN)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _EC.Sum() % _EC._pp % _EC._pu % _EC._uu).str()
        << flush << (format("  + EK0(FGC->BGP)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _E0.Sum() % _E0._pp % _E0._pu % _E0._uu).str() 
        << flush << (format("  + EDQ(FGC->MGN)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _EDQ.Sum() % _EDQ._pp % _EDQ._pu % _EDQ._uu).str()
        << flush << (format("  + EJ(shape-dep.) = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _EJ.Sum() % _EJ._pp % _EJ._pu % _EJ._uu).str()
        << flush << (format("  = -----------------------------------------------------------------------------------")).str()
        << flush << (format("  + SUM(E) (0,Q,J) = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _ET.Sum() % _ET._pp % _ET._pu % _ET._uu).str()
        << flush;
    
    double inner_epp = _polar_EPP;
    double inner_eppu = _polar_EF00+_polar_EF01+_polar_EF02+_polar_EF11+_polar_EF12 - _polar_EPP;
    double inner_ework = _polar_EM0+_polar_EM1;
    double inner = inner_epp+inner_eppu+inner_ework;
    
    LOG(logINFO,*_log)
        << (format("Interaction FGC <> FGC")).str()
        << flush << (format("  + EF [00 01 11]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % (_polar_EF00+_polar_EF01+_polar_EF11) % _polar_EF00 % _polar_EF01 % _polar_EF11).str()
        << flush << (format("  + EF [02 12 22]  = %1$+1.7e = *****ZERO*****  *****ZERO*****  *****ZERO***** eV") 
            % 0.0).str()
        << flush << (format("  + EM [0  1  2 ]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  *****ZERO***** eV") 
            % (_polar_EM0+_polar_EM1) % _polar_EM0 % _polar_EM1).str()
        << flush << (format("  o E  [PP PU UU]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % (_polar_EPP+_polar_EPU+_polar_EUU) % _polar_EPP % _polar_EPU % _polar_EUU).str()
        << flush << (format("  = -----------------------------------------------------------------------------------")).str()
        << flush << (format("  + SUM(E)         = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV")
            % inner % inner_epp % inner_eppu % inner_ework).str()
        << flush;
    
    _Estat = outer_epp + inner_epp;
    _Eindu = outer_eppu + inner_eppu + inner_ework;
    _Eppuu = _Estat + _Eindu;
    LOG(logINFO,*_log)
            
        << (format("Interaction FGC <> FGC(i) u ***(o)")).str()
        << flush << (format("  + Ei [pp+pu+iw]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV")
            % inner % inner_epp % inner_eppu % inner_ework).str()
        << flush << (format("  + Eo [pp+pu+iw]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  ************** eV")
            % outer % outer_epp % outer_eppu).str()
        << flush << (format("  = ===================================================================================")).str()
        << flush << (format("  + E  [stat+ind]  = %1$+1.7e = %2$+1.7e  %3$+1.7e eV")
            % _Eppuu % _Estat % _Eindu).str()
        << flush;
    
    
    _t_total = _t_coarsegrain+_t_fields+_t_induction+_t_energy;
    
    LOG(logDEBUG,*_log) << flush << (format("Timing (T = %1$1.2f min)") % (_t_total)) << flush;
    LOG(logDEBUG,*_log) << (format("  o Usage <Coarsegrain> = %1$2.2f%%") % (100*_t_coarsegrain/_t_total)) << flush;
    LOG(logDEBUG,*_log) << (format("  o Usage <Fields>      = %1$2.2f%%") % (100*_t_fields/_t_total)) << flush;
    LOG(logDEBUG,*_log) << (format("  o Usage <Induction>   = %1$2.2f%%") % (100*_t_induction/_t_total)) << flush;
    LOG(logDEBUG,*_log) << (format("  o Usage <Energy>      = %1$2.2f%%") % (100*_t_energy/_t_total)) << flush;    
    LOG(logDEBUG,*_log) << flush;
    
    for (vector<PolarSeg*>::iterator sit1 = _fg_C.begin(); 
        sit1 != _fg_C.end(); ++sit1) {
        (*sit1)->ClearPolarNbs();
    }    
    return;
}


void Ewald3DnD::EvaluateFields() {
    
    // RESET PERMANENT FIELDS
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }

    // REAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    Field_ConvergeRealSpaceSum();    

    // RECIPROCAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    Field_ConvergeReciprocalSpaceSum();

    // SHAPE-CORRECTION (3D3D)/ K0-CORRECTION (3D2D)
    Field_CalculateShapeCorrection();

    // FOREGROUND CORRECTION (3D2D && 3D3D)
    Field_CalculateForegroundCorrection();
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    LOG(logDEBUG,*_log) << flush << "Foreground fields:" << flush;
    int fieldCount = 0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {        
        PolarSeg* pseg = *sit1;        
        for (pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
            vec fp = (*pit1)->getFieldP();
            LOG(logDEBUG,*_log)
               << (format("F = (%1$+1.7e %2$+1.7e %3$+1.7e) V/m") 
                    % (fp.getX()*EWD::int2V_m) 
                    % (fp.getY()*EWD::int2V_m) 
                    % (fp.getZ()*EWD::int2V_m)).str() << flush;
            fieldCount += 1;
            if (fieldCount > 10) {
                LOG(logDEBUG,*_log)
                    << "F = ... ... ..." << flush;
                break;
            }
        }
        if (fieldCount > 10) break;
    }
    
    return;
}


void Ewald3DnD::EvaluateInduction() {
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << format("Call inductor on FGC = QM0 u MM1") << flush;
    LOG(logDEBUG,*_log) << (format("  o |QM0|, |MM1|, |MM2|        %1$d %2$d %3$d") 
            % _ptop->QM0().size() % _ptop->MM1().size() % _ptop->MM2().size()).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Polarization cut-off:      ")).str() << _polar_cutoff << " nm " << flush;
    LOG(logDEBUG,*_log) << (format("  o With induction:            %1$s") % ((_polar_do_induce) ? "yes" : "no")) << flush;
    LOG(logDEBUG,*_log) << (format("  o Thole sharpness parameter: ")).str() << _polar_aDamp << flush;
    LOG(logDEBUG,*_log) << (format("  o SOR mixing factor:         ")).str() << _polar_wSOR_N << " (N) " << _polar_wSOR_C << " (C) "  << flush;
    LOG(logDEBUG,*_log) << (format("  o Iterations (max):          512")).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Tolerance (rms, e*nm):     0.001")).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Induce within QM0:         yes")).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Subthreads:                single")).str() << flush;
    
    // Forge XJob object to comply with XInductor interface
    bool polar_has_permanent_fields = true;
    XJob polar_xjob = XJob(_ptop, polar_has_permanent_fields);
    
    // INITIALIZE XINDUCTOR
    bool    polar_induce_intra_pair = true;
    int     polar_subthreads = 1;
    double  polar_epstol = 0.001;
    int     polar_maxIter = 512;
    bool    polar_maverick = _log->isMaverick(); // TODO Extract from _log
    
    XInductor polar_xind = XInductor(_polar_do_induce, 
                                     polar_induce_intra_pair, 
                                     polar_subthreads,
                                     _polar_wSOR_N,
                                     _polar_wSOR_C,
                                     polar_epstol,
                                     polar_maxIter,
                                     _polar_aDamp,
                                     polar_maverick,
                                     _top);
    polar_xind.setLog(_log);
    polar_xind.Evaluate(&polar_xjob);
    
    // SAVE CONVERGENCE
    _polar_converged = polar_xind.hasConverged();
    
    // SAVE RESULTS
    _polar_ETT = polar_xjob.getETOT();
    _polar_EPP = polar_xjob.getEPP();
    _polar_EPU = polar_xjob.getEPU();
    _polar_EUU = polar_xjob.getEUU();
    
    _polar_EF00 = polar_xjob.getEF00();
    _polar_EF01 = polar_xjob.getEF01();
    _polar_EF02 = polar_xjob.getEF02();
    _polar_EF11 = polar_xjob.getEF11();
    _polar_EF12 = polar_xjob.getEF12();
    _polar_EM0 = polar_xjob.getEM0();
    _polar_EM1 = polar_xjob.getEM1();
    _polar_EM2 = polar_xjob.getEM2();
    return;
}
        
        
void Ewald3DnD::EvaluateEnergy() {
    
    // REAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    EWD::triple<> EPP_fgC_mgN = ConvergeRealSpaceSum();    
    
    // RECIPROCAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    EWD::triple<> EKK_fgC_bgP = ConvergeReciprocalSpaceSum();       
    
    // K=0 TERM (FOR 3D2D)
    EWD::triple<> EK0_fgC_bgP = CalculateK0Correction();
    
    // SHAPE-CORRECTION (FOR 3D3D)
    EWD::triple<> EJ_fgC_bgP = CalculateShapeCorrection();    
    
    // REAL-SPACE HIGHER-RANK CORRECTION (3D2D && 3D3D)
    EWD::triple<> EDQ_fgC_mgN = CalculateHigherRankCorrection();
    
    // FOREGROUND CORRECTION (3D2D && 3D3D)
    EWD::triple<> EPP_fgC_fgN = CalculateForegroundCorrection();    
    
    double int2eV = EWD::int2eV;
    _ER  = EPP_fgC_mgN * int2eV;
    _EK  = EKK_fgC_bgP * int2eV;
    _E0  = EK0_fgC_bgP * int2eV;
    _EJ  = EJ_fgC_bgP  * int2eV;
    _EDQ = EDQ_fgC_mgN * int2eV;
    _EC  = EPP_fgC_fgN * int2eV;
    _ET  = _ER + _EK + _E0 + _EJ + _EDQ - _EC;
    
    return;
}


void Ewald3DnD::EvaluatePoisson() {
    
    POI::PoissonGrid poisson_grid(_top, _fg_C, _bg_P, _log);
}


EWD::triple<> Ewald3DnD::ConvergeRealSpaceSum() {
    
    LOG(logDEBUG,*_log) << flush;

    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    // REAL-SPACE CONVERGENCE   
    double dR = 0.1; // radius increment [nm]
    _converged_R = false;
    double prev_ER = 0.0;
    double this_ER = 0.0;
    for (int i = 0; i < 1000; ++i) {
        double Rc = _R_co + (i-1)*dR;
        // Set-up midground
        this->SetupMidground(Rc);
        // Calculate interaction energy
        this_ER = 0.;
        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
            for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        this_ER += _actor.E_QQ_ERFC(*(*pit1), *(*pit2), _alpha);
                    }
                }
            }
        }
        double dER_rms = sqrt((this_ER-prev_ER)*(this_ER-prev_ER))*EWD::int2eV;
        LOG(logDEBUG,*_log)
            << (format("Rc = %1$+1.7f   |MGN| = %3$5d nm   ER = %2$+1.7f eV   dER(rms) = %4$+1.7f") 
            % Rc % (this_ER*EWD::int2eV) % _mg_N.size() % dER_rms).str() << flush;
        if (i > 0 && dER_rms < _crit_dE) {
            _converged_R = true;
            LOG(logDEBUG,*_log)  
                << (format(":::: Converged to precision as of Rc = %1$+1.3f nm") 
                % Rc ) << flush;
            break;
        }
        prev_ER = this_ER;
    }
    return EWD::triple<>(this_ER,0,0);
}


EWD::triple<> Ewald3DnD::CalculateForegroundCorrection() {
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    double EPP_fgC_fgN = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_fgN += _actor.E_QQ_ERF(*(*pit1), *(*pit2), _alpha);
                }
            }
        }
    }
    return EWD::triple<>(EPP_fgC_fgN,0,0);
}


EWD::triple<> Ewald3DnD::CalculateHigherRankCorrection() {
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    double EDQ_fgC_mgN = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    _actor.BiasStat(*(*pit1), *(*pit2));
                    EDQ_fgC_mgN += _actor.E_Q0_DQ(*(*pit1), *(*pit2));
                }
            }
        }
    }
    return EWD::triple<>(EDQ_fgC_mgN,0,0);
}


string Ewald3DnD::GenerateErrorString() {
    string rstr;
    rstr += (format("Converged R-sum = %1$s, converged K-sum = %2$s, ")
        % ((_converged_R) ? "true" : "false")
        % ((_converged_K) ? "true" : "false")).str();
    rstr += (format("converged induction = %1$s")
        % ((_polar_converged) ? "true" : "false")).str();
    return rstr;
}


Property Ewald3DnD::GenerateOutputString() {
    
    Property prop;
    Property &out = prop.add("output","");
    Property *next = NULL;    
    
    next = &out.add("summary", "");
    next->add("type", _jobType);
    next->add("xyz", (format("%1$+1.7f %2$+1.7f %3$+1.7f") 
        % _center.getX() % _center.getY() % _center.getZ()).str())
        .setAttribute("unit","nm");
    next->add("total", (format("%1$+1.7f") 
        % _Eppuu).str())
        .setAttribute("unit","eV");
    next->add("estat", (format("%1$+1.7f") 
        % _Estat).str())
        .setAttribute("unit","eV");
    next->add("eindu", (format("%1$+1.7f") 
        % _Eindu).str())
        .setAttribute("unit","eV");
    
//    next = &out.add("splitting", "");
//    next->add("R-term", (format("%1$+1.7f") % _ER.Sum()).str());
//    next->add("K-term", (format("%1$+1.7f") % _EK.Sum()).str());
//    next->add("O-term", (format("%1$+1.7f") % _E0.Sum()).str());
//    next->add("J-term", (format("%1$+1.7f") % _EJ.Sum()).str());
//    next->add("C-term", (format("%1$+1.7f") % _EC.Sum()).str());
//    next->add("Q-term", (format("%1$+1.7f") % _EDQ.Sum()).str());
    
    next = &out.add("terms_i", "");
    next->add("F-00-01-11", (format("%1$+1.5e %2$+1.5e %3$+1.5e") % _polar_EF00 % _polar_EF01 % _polar_EF11).str());
    next->add("M-00-11---", (format("%1$+1.5e %2$+1.5e") % _polar_EM0 % _polar_EM1).str());
    next->add("E-PP-PU-UU", (format("%1$+1.5e %2$+1.5e %3$+1.5e") % _polar_EPP % _polar_EPU % _polar_EUU).str());
    
    next = &out.add("terms_o", "");
    next->add("R-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _ER.Sum() % _ER._pp % _ER._pu % _ER._uu).str());
    next->add("K-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _EK.Sum() % _EK._pp % _EK._pu % _EK._uu).str());
    next->add("O-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _E0.Sum() % _E0._pp % _E0._pu % _E0._uu).str());
    next->add("J-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _EJ.Sum() % _EJ._pp % _EJ._pu % _EJ._uu).str());
    next->add("C-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _EC.Sum() % _EC._pp % _EC._pu % _EC._uu).str());
    next->add("Q-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _EDQ.Sum() % _EDQ._pp % _EDQ._pu % _EDQ._uu).str());
    
    next = &out.add("shells", "");
    next->add("FGC", (format("%1$d") % _fg_C.size()).str());
    next->add("FGN", (format("%1$d") % _fg_N.size()).str());
    next->add("MGN", (format("%1$d") % _mg_N.size()).str());
    next->add("BGN", (format("%1$d") % _bg_N.size()).str());
    next->add("BGP", (format("%1$d") % _bg_P.size()).str());
    next->add("QM0", (format("%1$d") % _polar_qm0.size()).str());
    next->add("MM1", (format("%1$d") % _polar_mm1.size()).str());
    next->add("MM2", (format("%1$d") % _polar_mm2.size()).str());
    
    next = &out.add("timing", "");
    next->add("t_total", (format("%1$1.2f") % _t_total).str())
        .setAttribute("unit","min");
    next->add("t_wload", (format("%1$1.2f %2$1.2f %3$1.2f %4$1.2f")
        % (_t_coarsegrain) % (_t_fields)
        % (_t_induction)   % (_t_energy)).str())
        .setAttribute("unit","min");
        
    return prop;
}
    
    
    
    
    
    
}}
