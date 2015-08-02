#ifndef VOTCA_CTP_LOG2MPS_H
#define VOTCA_CTP_LOG2MPS_H


#include <votca/ctp/qmtool.h>
#include <votca/ctp/topology.h>
#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/qmmachine.h>


namespace votca { namespace ctp {

class Log2Mps : public QMTool
{
public:

    Log2Mps() { };
   ~Log2Mps() { };

    string Identify() { return "log2mps"; }

    void   Initialize(Property *options);
    bool   Evaluate();

    
private:

    string _package;
    string _logfile;
    string _mpsfile;
    
    
    //
    // start REDUCED md2qm
    // 
    bool   _mps_reduced;
    
    // containers    
    vector<int>    _mps_red_inds;
    vector<string> _mps_red_ats;
    vector<double> _mps_red_qs;
    
    // internal methods
    void Error_missing_opt(string opt_name){
        // flip an error on missing option
        cout << "\n\n\t" << "Error:" << endl;
        cout << "\t"     << "I can't find option <" << opt_name << ">, I need it" << endl;
        cout << "\n"     << endl;
        throw std::runtime_error("yes, it's an error");
    }
    
    double charge_sum(vector<QMAtom*>& qmats){
        // calc total charge 
        double Qtot = 0;
        for (int i=0; i < qmats.size(); ++i) Qtot += qmats[i]->charge;
        return Qtot;
    }
    
    bool superposition(vector<QMAtom*>& qmatoms,
                       vector<int>&     kill_list ){
        // brings qm and mm parts together
        // enjoy reading this
        sort(kill_list.begin(), kill_list.end());
        reverse(kill_list.begin(), kill_list.end());
        
        int finger(-1);
        vector<int>::iterator it_kl = kill_list.begin();
        for ( it_kl ; it_kl != kill_list.end(); ++it_kl ){
            // pick atom from the list, get it's charge, 
            // delete instance, delete pointer
            finger = (*it_kl)-1;
            double dq = (*( qmatoms.begin()+finger ))->charge;
            delete *( qmatoms.begin()+finger );
            qmatoms.erase( qmatoms.begin()+finger );
            // insert backwards missing atoms
            for (int i=_mps_red_ats.size()-1; i!=-1; --i)
                qmatoms.insert( qmatoms.begin()+finger,
                    new QMAtom( _mps_red_ats[i], -1, -1, -1, _mps_red_qs[i], false ) );
            // first one gets the qm charge of the erased atom
            (*( qmatoms.begin()+finger ))->charge += dq;

        }
        return true;
    }
    
};


void Log2Mps::Initialize(Property *opt) {
    
    QMPackageFactory::RegisterAll();
    
    string key = "options.log2mps";
    _package = opt->get(key+".package").as<string>();
    _logfile = opt->get(key+".logfile").as<string>();
    
    
    _mpsfile = (opt->exists(key+".mpsfile")) ? 
        opt->get(key+".mpsfile").as<string>() : "";
    if (_mpsfile == "") _mpsfile = _logfile.substr(0,_logfile.size()-4)+".mps";
    
    //
    // start REDUCED md2qm
    // 
    _mps_reduced = (opt->exists(key+".reduced"));
    if (_mps_reduced){
        string inds, atoms, charges;
        string key2 = key+".reduced";
        
        // read in opt strings
        if (opt->exists(key2+".inds"))    
            inds    = opt->get( key2+".inds"   ).as<string>();
        else 
            Error_missing_opt("inds");
           
        if (opt->exists(key2+".atoms"))   
            atoms   = opt->get( key2+".atoms"  ).as<string>();
        else 
            Error_missing_opt("atoms");

        if (opt->exists(key2+".charges")) 
            charges = opt->get( key2+".charges").as<string>();
        else 
            Error_missing_opt("charges");
        
        // convert strings to internal data types
        // todo | no convert handling
        // todo | equal qs and ats handling
        // todo | tot charge of md = 0 handling
        Tokenizer toker_charges(charges, " \n");
        toker_charges.ConvertToVector<double>(_mps_red_qs );
        
        Tokenizer toker_inds(   inds,    " \n");
        toker_inds.ConvertToVector<int>(     _mps_red_inds);
        
        Tokenizer toker_atoms(  atoms,   " \n");
        toker_atoms.ToVector(                _mps_red_ats );
        
    } // END of REDUCED
    
    cout << endl << "... ... " << _logfile << " => " << _mpsfile << flush;
}


bool Log2Mps::Evaluate() {
    
    // Logger (required for QM package, so we can just as well use it)
    Logger log;
    log.setPreface(logINFO, "\n... ...");
    log.setPreface(logDEBUG, "\n... ...");
    log.setReportLevel(logDEBUG);
    log.setMultithreading(true);  
    
    // Set-up QM package
    LOG(logINFO,log) << "Using package <" << _package << ">" << flush;
    QMPackage *qmpack = QMPackages().Create(_package);    
    qmpack->doGetCharges(true);
    qmpack->setLog(&log);
    qmpack->setRunDir(".");
    qmpack->setLogFileName(_logfile);
    
    // Create orbitals, fill with life & extract QM atoms
    Orbitals orbs;
    int cdx = qmpack->ParseLogFile(&orbs);
    if (!cdx) {
        cout << "\nERROR Parsing " << _logfile << "failed. Abort." << endl;
        throw std::runtime_error("(see above, parsing error)");
    }    
    vector<QMAtom*> &qmatoms = *orbs.getAtoms();
    vector<QMAtom*>::iterator it;
    
    if (qmatoms.size() < 1) {
        cout << "\nERROR No charges extracted from " << _logfile 
            << ". Abort.\n" << flush;
        throw std::runtime_error("(see above, input or parsing error)");
    }
        
    //
    // start REDUCED md2qm
    // 
    if (_mps_reduced){
        // superimpose qm and missing md
        superposition(qmatoms,_mps_red_inds);
    } // END of REDUCED
    
    // Sanity checks, total charge
    LOG(logINFO,log) << 
            qmatoms.size() << " QM atoms, total charge Q = " << 
                charge_sum(qmatoms) << flush;    
    
    // Convert to polar segment & write mps-file
    QMMInterface qmmface;
    PolarSeg *pseg = qmmface.Convert(qmatoms);
    
    string tag = "::LOG2MPS " 
        + (boost::format("(log-file='%1$s' : %2$d QM atoms)")
        % _logfile % qmatoms.size()).str();    
    pseg->WriteMPS(_mpsfile, tag);
    return true;
}



}}

#endif