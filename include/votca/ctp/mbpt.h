/*
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// UBLAS stops checking types and array bounds if this flag is defined
#define NDEBUG
#define BOOST_UBLAS_NDEBUG

#ifndef _VOTCA_CTP_MBPT_H
#define	_VOTCA_CTP_MBPT_H

#include <votca/ctp/segment.h>
#include <votca/ctp/orbitals.h>
#include <votca/ctp/aobasis.h>
#include <votca/ctp/aomatrix.h>
#include <votca/ctp/threecenters.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/parallelxjobcalc.h>
#include <unistd.h>

#include <fstream>
#include <sys/stat.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/linalg.h>

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
// #include <gsl/gsl_eigen.h>
// #include <gsl/gsl_linalg.h>
// #include <gsl/gsl_cblas.h>

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
/**
* \brief GWBSE implementation
*
* Requires a first-principles package, i.e. GAUSSIAN installation
*
* Callname: gwbse
*/

class MBPT 
{
public:

    MBPT() {};
   ~MBPT() {};

  /*  string  Identify() { return "gwbse"; }
    void    Initialize( Property *options);
    void    ParseOrbitalsXML(Topology *top, Property *options);
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *thread);
 */
    void    CleanUp();

    // int getMlower(){ return mmin -1; };
    // int getMupper(){ return mmax -1; };
    
    void setLogger( Logger* pLog ) { _pLog = pLog; }
    
    bool Evaluate(   Orbitals* _orbitals );

    // interfaces for options getting/setting
    bool get_do_qp_diag(){ return _do_qp_diag ;}
    void set_do_qp_diag( bool inp ){ _do_qp_diag = inp;}

    bool get_do_bse_singlets(){ return _do_bse_singlets ;}
    void set_do_bse_singlets( bool inp ){ _do_bse_singlets = inp;}

    bool get_do_bse_triplets(){ return _do_bse_triplets ;}
    void set_do_bse_triplets( bool inp ){ _do_bse_triplets = inp;}

    bool get_store_qp_pert(){ return _store_qp_pert ;}
    void set_store_qp_pert( bool inp ){ _store_qp_pert = inp;}
    
    bool get_store_qp_diag(){ return _store_qp_diag ;}
    void set_store_qp_diag( bool inp ){ _store_qp_diag = inp;}
    
    bool get_store_bse_singlets(){ return _store_bse_singlets ;}
    void set_store_bse_singlets( bool inp ){ _store_bse_singlets = inp;}
    
    bool get_store_bse_triplets(){ return _store_bse_triplets ;}
    void set_store_bse_triplets( bool inp ){ _store_bse_triplets = inp;}
    
    string get_ranges(){ return _ranges ;}
    void set_ranges( string inp ){ _ranges = inp;}
    
    double get_rpamaxfactor() {return _rpamaxfactor ;}
    void set_rpamaxfactor( double inp ) { _rpamaxfactor = inp; }
  
    double get_qpmaxfactor() {return _qpmaxfactor ;}
    void set_qpmaxfactor( double inp ) { _qpmaxfactor = inp; }

    double get_qpminfactor() {return _qpminfactor ;}
    void set_qpminfactor( double inp ) { _qpminfactor = inp; }

    double get_bsemaxfactor() {return _bsemaxfactor ;}
    void set_bsemaxfactor( double inp ) { _bsemaxfactor = inp; }

    double get_bseminfactor() {return _bseminfactor ;}
    void set_bseminfactor( double inp ) { _bseminfactor = inp; }
    
    unsigned int get_bse_vmin() {return _bse_vmin ;}
    void set_bse_vmin( unsigned int inp ) { _bse_vmin = inp; }

    unsigned int get_bse_cmax() {return _bse_cmax ;}
    void set_bse_cmax( unsigned int inp ) { _bse_cmax = inp; }

    unsigned int get_rpamax() {return _rpamax ;}
    void set_rpamax( unsigned int inp ) { _rpamax = inp; }

    unsigned int get_qpmax() {return _qpmax ;}
    void set_qpmax( unsigned int inp ) { _qpmax = inp; }

    unsigned int get_qpmin() {return _qpmin ;}
    void set_qpmin( unsigned int inp ) { _qpmin = inp; }
    
    int get_bse_nmax(){return _bse_nmax;}
    void set_bse_nmax( int inp){ _bse_nmax = inp;}
    
    string get_gwbasis_name(){return _gwbasis_name;}
    void set_gwbasis_name(string inp){ _gwbasis_name = inp;}
    
    string get_dftbasis_name(){return _dftbasis_name;}
    void set_dftbasis_name(string inp){ _dftbasis_name = inp;}

    double get_shift() {return _shift ;}
    void set_shift( double inp ) { _shift = inp; }
    
    private:

    Logger *_pLog;
    

    
    //bool   _maverick;
    
    // program tasks
    bool                                _do_qp_diag;
    bool                                _do_bse_singlets;
    bool                                _do_bse_triplets;
    
    // storage tasks
    bool                                _store_qp_pert;
    bool                                _store_qp_diag;
    bool                                _store_bse_singlets;
    bool                                _store_bse_triplets;
    
    
    
    string _outParent;
    string _outMonDir;
    
    string _package;
    Property _package_options;   
    
    string _gwpackage;
    Property _gwpackage_options; 
    
    // basis sets
    string                              _gwbasis_name;
    string                              _dftbasis_name;

    string                              _ranges;          // range types
    unsigned int                        _homo;            // HOMO index
    unsigned int                        _rpamin;
    unsigned int                        _rpamax;
    double                              _rpamaxfactor;    // RPA level range
    unsigned int                        _qpmin;
    unsigned int                        _qpmax;
    unsigned int                        _qptotal;
    double                              _qpminfactor;
    double                              _qpmaxfactor;     // QP level range
    double                              _bseminfactor;
    double                              _bsemaxfactor;
    unsigned int                        _bse_vmin;
    unsigned int                        _bse_vmax;
    unsigned int                        _bse_cmin;
    unsigned int                        _bse_cmax;
    unsigned int                        _bse_size;
    unsigned int                        _bse_vtotal;
    unsigned int                        _bse_ctotal;
    int                                 _bse_nmax;
         
    double                              _shift;  // pre-shift of DFT energies

    
    // RPA related variables and functions
    // container for the epsilon matrix
    std::vector< ub::matrix<double> > _epsilon;
    // container for frequencies in screening (index 0: real part, index 1: imaginary part)
    ub::matrix<double> _screening_freq;
    void symmetrize_threecenters(TCMatrix& _Mmn, ub::matrix<double>& _coulomb);
    void RPA_calculate_epsilon( TCMatrix& _Mmn_RPA , ub::matrix<double> _screening_freq , double _shift , ub::vector<double>& _dft_energies  );
    void RPA_prepare_threecenters( TCMatrix& _Mmn_RPA, TCMatrix& _Mmn_full, AOBasis& gwbasis, AOMatrix& gwoverlap, AOMatrix& gwoverlap_inverse     );

    
    // PPM related variables and functions
    ub::matrix<double> _ppm_phi;
    ub::vector<double> _ppm_freq;
    ub::vector<double> _ppm_weight;
    
    void PPM_construct_parameters( ub::matrix<double>& _overlap_cholesky_inverse   );
    
    // Sigma related variables and functions
    ub::matrix<double> _sigma_x; // exchange term
    ub::matrix<double> _sigma_c; // correlation term
    
    void sigma_prepare_threecenters( TCMatrix& _Mmn );
    void sigma_x_setup(const TCMatrix& _Mmn );
    void sigma_c_setup(const TCMatrix& _Mmn , const ub::vector<double>& _edft );
    
    // QP variables and functions
    ub::vector<double> _qp_energies;
    ub::matrix<double> _vxc;
    ub::vector<double> _qp_diag_energies;     // those should be directly stored in 
    ub::matrix<double> _qp_diag_coefficients; // orbitals object, once the interface is set
    void FullQPHamiltonian();
    
    // BSE variables and functions
    ub::matrix<double> _eh_x;
    ub::matrix<double> _eh_d;
    ub::matrix<double> _eh_qp;
    ub::vector<double> _bse_singlet_energies;
    ub::matrix<double> _bse_singlet_coefficients;
    ub::vector<double> _bse_triplet_energies;
    ub::matrix<double> _bse_triplet_coefficients;
    
    std::vector< ub::matrix<double> > _interlevel_dipoles;
    std::vector< ub::matrix<double> > _interlevel_dipoles_electrical;
    void BSE_x_setup(const TCMatrix& _Mmn );
    void BSE_d_setup(const TCMatrix& _Mmn );
    void BSE_qp_setup( );
    void BSE_solve_triplets();
    void BSE_solve_singlets();
    
};


}}

#endif	/* _VOTCA_CTP_MBPT_H */