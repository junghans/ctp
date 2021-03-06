#ifndef VOTCA_CTP_EWDBGPOL_H
#define VOTCA_CTP_EWDBGPOL_H

#include <votca/ctp/qmthread.h>
#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/polarbackground.h>

namespace votca { 
namespace ctp {
  
class EwaldBgPolarizer : public QMCalculator
{
public:

    EwaldBgPolarizer() {};
   ~EwaldBgPolarizer() {};

    string Identify() { return "ewdbgpol"; }
    void Initialize(Property *options);
    bool EvaluateFrame(Topology *top);

private:

    Property                       *_options;
    string                         _mps_table;
    string                         _xml_file;
    XMpsMap                        _mps_mapper;
    bool                           _pdb_check;

};



void EwaldBgPolarizer::Initialize(Property *opt) {
    
    _options = opt;
    _maverick = (_nThreads == 1) ? true : false;
    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;    

    string key = "options.ewdbgpol";
        if (opt->exists(key+".multipoles")) {
            _xml_file = opt->get(key+".multipoles").as< string >();
        }
        else {
            cout << endl;
            throw std::runtime_error("No multipole mapping file provided");
        }
    key = "options.ewdbgpol.control";
        if (opt->exists(key+".mps_table")) {
            _mps_table = opt->get(key+".mps_table").as<string>();
        }
        else {
            _mps_table = opt->get(key+".emp_file").as<string>();
        }
        if (opt->exists(key+".pdb_check")) {
            _pdb_check = opt->get(key+".pdb_check").as<bool>();
        }
        else { _pdb_check = false; }

    return;
}




bool EwaldBgPolarizer::EvaluateFrame(Topology *top) {

    QMThread master;
    master.getLogger()->setReportLevel(logDEBUG);
    master.getLogger()->setMultithreading(true);
    master.getLogger()->setPreface(logINFO,    "\nMST INF");
    master.getLogger()->setPreface(logERROR,   "\nMST ERR");
    master.getLogger()->setPreface(logWARNING, "\nMST WAR");
    master.getLogger()->setPreface(logDEBUG,   "\nMST DBG");
    Logger &log = *master.getLogger();
    
    
    // GENERATE BACKGROUND (= periodic bg, with empty foreground)
    cout << endl << "... ... Initialize MPS-mapper: " << flush;
    PolarTop ptop(top);
    _mps_mapper.GenerateMap(_xml_file, _mps_table, top);
    _mps_mapper.Gen_BGN(top, &ptop, &master);
    if (_pdb_check) ptop.PrintPDB("ewdbgpol.ptop.pdb");
    
    // POLARIZE SYSTEM
    EWD::PolarBackground pbg(top, &ptop, _options, &log);
    pbg.Polarize(_nThreads);
    
    // SAVE POLARIZATION STATE
    LOG(logINFO,log) << "Save polarization state" << flush;
    ptop.SaveToDrive("bgp_main.ptop");
    ptop.PrintPDB("bgp_main.pdb");
    
//    // LOAD POLARIZATION STATE
//    LOG(logINFO,log) << "Load polarization state" << flush;
//    PolarTop ptop2;
//    ptop2.LoadFromDrive("bgp_main.ptop");
//    ptop2.PrintPDB("bgp_check.pdb");
//    ptop2.SaveToDrive("bgp_check.ptop");
    
    return true;
}
    
}}

#endif