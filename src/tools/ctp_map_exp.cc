#include <iostream>
#include <fstream>
#include <stdexcept>

#include <votca/ctp/qmtopology.h>
#include <votca/csg/csgapplication.h>
#include "votca/tools/application.h"
#include <votca/csg/trajectorywriter.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/csg/topologyreader.h>
 
using namespace std;
using namespace votca::csg;
using namespace votca::ctp;
using namespace votca::tools;



class CtpMapExp : public Application
{

public:
    string ProgramName()  { return "ctp_map_exp"; }
    void   HelpText(ostream &out) {out << "Generates QM|MD topology" << endl;}

    void Initialize();
    void Run();
    bool EvaluateOptions();
    void BeginEvaluate() { ; }
    bool DoTrajectory() { return 0; }
    bool DoMapping() { return 0; }


protected:
    Property       _options;
    Topology       _mdtopol;
    Topology       _qmtopol;

};

namespace propt = boost::program_options;

void CtpMapExp::Initialize() {

    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();

    AddProgramOptions() ("top", propt::value<string> (),
                         "  Atomistic topology file ");
    AddProgramOptions() ("trj", propt::value<string> (),
                         "  Atomistic trajetory file ");
    AddProgramOptions() ("cg",  propt::value<string> (),
                         "  Coarse-Graining definitions ");
    AddProgramOptions() ("segments,s", propt::value<string> (),
                         "  Conjugated-segment definitions ");
    AddProgramOptions() ("file,f", propt::value<string> (),
                         "  SQLite state file ");
}

bool CtpMapExp::EvaluateOptions() {

    CheckRequired("top", "Missing topology file");
    CheckRequired("cg", "Missing CG definition file");
    CheckRequired("segments");
    CheckRequired("file");

    return 1;
}

void CtpMapExp::Run() {

    string topfile = _op_vm["top"].as<string> ();
    TopologyReader *reader;
    reader = TopReaderFactory().Create(topfile);

    if (reader == NULL) {
        throw runtime_error( string("Input format not supporeted: ")
                           + _op_vm["top"].as<string> () );
    }

    reader->ReadTopology(topfile, this->_mdtopol);
    cout << "MD Topology from " << topfile << ": Found "
         << _mdtopol.BeadCount() << " atoms in "
         << _mdtopol.MoleculeCount() << " molecules. "
         << endl;
    
}


int main(int argc, char** argv)
{
    CtpMapExp ctpmap;
    return ctpmap.Exec(argc, argv);
 }
