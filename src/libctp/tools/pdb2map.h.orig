#ifndef _PDB2Map_H
#define _PDB2Map_H


#include <votca/ctp/topology.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/logger.h>
#include <boost/algorithm/string.hpp>
#include <votca/tools/vec.h>
#include <boost/format.hpp>

#include <votca/ctp/trajiofactory.h>

namespace votca { namespace ctp {
    namespace ba = boost::algorithm;
//    using namespace std;
    
class PDB2Map : public QMTool
{
public:

    PDB2Map() { };
   ~PDB2Map() { };
   
    string Identify() { return "pdb2map"; }
    
// reads options    
    void   Initialize(Property *options);
    
// makes xml
    bool   Evaluate();
    
// main functions
    void topMdQm2xml();   // generate XML from two above
    void adaptQM2MD();    // make QM and MD compatible
    void isConvertable(); // can PDB be used as QM ?
// secondary functions    
    void fileExtension(string &, const string &); // return file extension
    void setPeriodicTable(); // fill periodic table -> elemNumber:mass
    
// nice error function    
    void error1(string line){cout<<"\n\n"; throw runtime_error(line+'\n');}
    
private:

    string      _input_md_file;
    string      _input_qm_file;
    string      _output_file;
    
    bool        _conversion_possible;
    
    Topology    _MDtop;
    Topology    _QMtop;

   // element:mass map
    map <string,int> el2mass;

};

void PDB2Map::Initialize(Property* options) 
{   
// register trajIOs
    trajIOFactory::RegisterAll();
    
// can convert from MD to QM, if QM is not found
    _conversion_possible = false;
    
// read options    
    string key = "options.pdb2map.";
    
// md file
    if ( options->exists(key+"md") ){
        _input_md_file      = options->get(key+"md").as<string> ();
        cout << endl << "... ... MD input: \t" << _input_md_file;
    }
    else{
        error1("... ... <MD> option error.");
    }
    
// qm file
    if ( options->exists(key+"qm") ){
        _input_qm_file      = options->get(key+"qm").as<string> ();
        cout << endl << "... ... QM input: \t" << _input_qm_file;
    }
    else{
        string extension;
        fileExtension(extension,_input_md_file);
        
        if (extension == "pdb"){
         _input_qm_file = _input_md_file;
         _conversion_possible = true;
         cout << endl << "... ... QM input is not given. "
              << endl << "... ... Try conversion from MD."
              << endl << "... ... QM input (same as MD): \t" << _input_qm_file;
        }
        else{
            error1("... ... <QM> option error. Can convert only PDB.");
        }
    }

// output xml file
    if ( options->exists(key+"output") ){
        _output_file      = options->get(key+"output").as<string> ();
        cout << endl << "... ... Output: \t" << _output_file;
    }
    else{
        _output_file = "system_output.xml";
        cout << endl << "... ... Output by default: \t" << _output_file;
    }   
}

bool PDB2Map::Evaluate() {
    
    trajIO * trajPtr = 0;
    string extension;

    cout << endl << "... ... MD topology";
    
    fileExtension(extension,_input_md_file);
    trajPtr = trajIOs().Create(extension);
    trajPtr->read(_input_md_file,&_MDtop);
    
    cout << endl << "... ... QM topology";
    
    fileExtension(extension,_input_qm_file);
    trajPtr = trajIOs().Create(extension);
    trajPtr->read(_input_qm_file,&_QMtop);
    
    // if "can convert", try to convert
    // else XYZ is supplied, adapt it to MD
    if (_conversion_possible) isConvertable();
    else adaptQM2MD();
    
    topMdQm2xml();
    
////    LOG( logINFO, _log ) << "Reading from: " << _input_file << flush;    
////    cout << _log;
}

void PDB2Map::adaptQM2MD()
{
    
    // check if atom number is the same in QM and MD
    int numMDatoms = _MDtop.getMolecule(1)->NumberOfAtoms();
    int numQMatoms = _QMtop.getMolecule(1)->NumberOfAtoms();
    
    bool _QM2MDcompatible = (numMDatoms == numQMatoms) ? true : false;
    
    // if so, proceed
    if (_QM2MDcompatible){
        Molecule * MDmolecule = _MDtop.getMolecule(1);
        Molecule * QMmolecule = _QMtop.getMolecule(1);
        
        vector<Segment*> MDsegments = MDmolecule->Segments();
        vector<Segment*> QMsegments = QMmolecule->Segments();
        
        vector<Segment*>::iterator MDSegIt;
        vector<Segment*>::iterator QMSegIt;
        
        for(MDSegIt = MDsegments.begin(), QMSegIt = QMsegments.begin();
                MDSegIt < MDsegments.end();
                MDSegIt++, QMSegIt++  )
        {
            Fragment * QMfragment = 0;
            
            vector<Atom*> MDSegAtoms = (*MDSegIt)->Atoms();
            vector<Atom*> QMSegAtoms = (*QMSegIt)->Atoms();
            
            vector<Atom*>::iterator MDSegAtIt;
            vector<Atom*>::iterator QMSegAtIt;
            
            int old_res_num(-1),new_res_num(-1);
            string res_name = "bad_wolf";
            
            for(MDSegAtIt = MDSegAtoms.begin(), QMSegAtIt = QMSegAtoms.begin();
                    MDSegAtIt < MDSegAtoms.end();
                    MDSegAtIt++, QMSegAtIt++)
            {
                new_res_num = (*MDSegAtIt)->getResnr();
                if ( new_res_num != old_res_num)
                {
                    old_res_num = new_res_num;
                    res_name = (*MDSegAtIt)->getResname();
                    
                    QMfragment = _QMtop.AddFragment(res_name);
                    QMmolecule->AddFragment(QMfragment);
                    (*QMSegIt)->AddFragment(QMfragment);
                    
                    QMfragment->setTopology(&_QMtop);
                    QMfragment->setMolecule(QMmolecule);
                    QMfragment->setSegment(*QMSegIt);
                    
                    (*QMSegAtIt)->setFragment(QMfragment);
                    QMfragment->AddAtom(*QMSegAtIt);                    
                }
                else
                {
                    (*QMSegAtIt)->setFragment(QMfragment);
                    QMfragment->AddAtom(*QMSegAtIt);    
                }
            }
        }
    }
    else{
        error1("... ... Number of MD atoms is different from QM.\n"
               "... ... If it's the case of reduced molecule, "
                        " I need a map.\n"
               " ... ... Tags: map\n"
               " ... ... NOT IMPLEMENTED\n");
    }
} /* END adaptQM2MD */

void PDB2Map::topMdQm2xml()
{
    cout << endl << "... ... (A)merging XML from MD and QM topologies.";

    Molecule * MDmolecule = _MDtop.getMolecule(1);
    Molecule * QMmolecule = _QMtop.getMolecule(1);
    
    //
    // xml stuff
    //
    
    Property record;
    Property *ptopology_p = &record.add("topology","");
    Property *pmolecules_p = &ptopology_p->add("molecules","");
    Property *pmolecule_p = &pmolecules_p->add("molecule","");
    pmolecule_p->add("name","random_molecule_name");
    pmolecule_p->add("mdname","name_from_topol.top");
    Property *psegments_p = &pmolecule_p->add("segments","");
    Property *psegment_p = &psegments_p->add("segment","");
    psegment_p->add("name","random_segment_name");
    
    // qc data
    psegment_p->add("qmcoords","QC_FILES/your_file_with.xyz");
    psegment_p->add("orbitals","QC_FILES/your_file_with.fort7");
    psegment_p->add("basisset","INDO");
    
    // hole
    psegment_p->add("torbital_h","number_of_electrons,150");
    psegment_p->add("U_cC_nN_h","0.0000");
    psegment_p->add("U_nC_nN_h","0.0000");
    psegment_p->add("U_cN_cC_h","0.0000");
    
    // electron
    psegment_p->add("torbital_e","number_of_electrons+1,151");
    psegment_p->add("U_cC_nN_e","0.0000");
    psegment_p->add("U_nC_nN_e","0.0000");
    psegment_p->add("U_cN_cC_e","0.0000");

    // mps entry
    psegment_p->add("multipoles_n","MP_FILES/file_with.mps");
    psegment_p->add("multipoles_h","MP_FILES/file_with.mps");
    psegment_p->add("multipoles_e","MP_FILES/file_with.mps");
    
    psegment_p->add("map2md","0");
    
    // main body 
    Property *pfragments_p = &psegment_p->add("fragments","");
                                        
    vector < Segment * > allMdSegments = MDmolecule->Segments();
    vector < Segment * > allQmSegments = QMmolecule->Segments();
  
    vector < Segment * >::iterator segMdIt;
    vector < Segment * >::iterator segQmIt;
                                        
    for ( segMdIt = allMdSegments.begin(), 
                segQmIt = allQmSegments.begin();
            
          (allMdSegments.size() > allQmSegments.size()) ? 
              segMdIt < allMdSegments.end() :
              segQmIt < allQmSegments.end();
            
          segMdIt++, segQmIt++)
    {
        
        vector < Fragment * > allMdFragments = (*segMdIt)->Fragments();
        vector < Fragment * > allQmFragments = (*segQmIt)->Fragments();
                                        
        vector < Fragment * >::iterator fragMdIt;
        vector < Fragment * >::iterator fragQmIt;
                                        
        for ( fragMdIt = allMdFragments.begin() , 
                fragQmIt = allQmFragments.begin();
                
              (allMdFragments.size() > allQmFragments.size()) ?
                  fragMdIt < allMdFragments.end() :
                  fragQmIt < allQmFragments.end();
                
              fragMdIt++,fragQmIt++ )
        {
            string mapName;            
            stringstream mapMdAtoms;
            stringstream mapQmAtoms;
            stringstream mapMpoles;
            stringstream mapWeight;
            stringstream mapFrame;

            mapName      = (*fragMdIt)->getName() ;
            
            int localCounter = 0;
            vector < Atom * > allMdAtoms = (*fragMdIt)->Atoms();
            vector < Atom * > allQmAtoms = (*fragQmIt)->Atoms();

            vector < Atom * >::iterator atomMdIt;
            vector < Atom * >::iterator atomQmIt;
            for ( atomMdIt = allMdAtoms.begin(),
                    atomQmIt = allQmAtoms.begin();
                    
                  (allMdAtoms.size() > allQmAtoms.size()) ? 
                      atomMdIt < allMdAtoms.end() :
                      atomQmIt < allQmAtoms.end();
                  
                  atomMdIt++, atomQmIt++ )
            {
                if (atomMdIt < allMdAtoms.end())
                {
                        mapMdAtoms << boost::format("%=13s") % 
                                (boost::format("%s:%s:%s") 
                                   % (*atomMdIt)->getResnr()
                                   % (*atomMdIt)->getResname()
                                   % (*atomMdIt)->getName()).str() ;
                        
                                if (el2mass.find((*atomQmIt)->getElement()) 
                                        != el2mass.end())
                                {
                        mapWeight << boost::format("%=13i")
                                           % el2mass[(*atomQmIt)->getElement()];
                                }
                                else
                                {
                        mapWeight << boost::format("%=13i") % " ";        
                                }
                }
                else
                {
                        mapMdAtoms << boost::format("%=13s") %
                                (boost::format("%s:%s:%s") 
                                        % " " % " " % " " ).str();
                        
                        mapWeight << boost::format("%=13i") % " " ;
                }
                
                if (atomQmIt < allQmAtoms.end())
                {
                        mapQmAtoms << boost::format("%=13s") %
                                (boost::format("%1%:%2% ") 
                                   % (*atomQmIt)->getId()
                                   % (*atomQmIt)->getElement()).str() ;
                }
                else
                {
                        mapQmAtoms << boost::format("%=13s") %
                                (boost::format("%1%:%2% ") 
                                        % " " % " " ).str();
                }
                
                if (localCounter < 3 && localCounter < allMdAtoms.size())
                {
                        mapFrame << boost::format("%=5i")
                                   % (*atomMdIt)->getId();
                }
                localCounter++;
            }
            
            mapMpoles << " " << mapQmAtoms.str();
            
            Property *pfragment_p  = &pfragments_p->add("fragment","");
            pfragment_p->add("name", mapName);
            pfragment_p->add("mdatoms", mapMdAtoms.str());
            pfragment_p->add("qmatoms", mapQmAtoms.str());
            pfragment_p->add("mpoles",  mapMpoles.str());
            pfragment_p->add("weights", mapWeight.str());
            pfragment_p->add("localframe", mapFrame.str());

         }
    }
    
//    cout << endl << setlevel(1) << XML << record;
        
    ofstream outfile( _output_file.c_str() );
    
    PropertyIOManipulator XML(PropertyIOManipulator::XML,1,""); 
    outfile << XML << record;
    outfile.close();
    
    cout << endl << "... ... Topology written to: " + _output_file;
    cout << endl << "... ... PDB2MAP is done. Victory.";
    
    return;
} /* END topMdQm2xml */

void PDB2Map::fileExtension(string & extension, const string & filename)
{
    Tokenizer tokenized_line( filename, ".");
    vector<string> vector_line;
    tokenized_line.ToVector(vector_line);
    extension = vector_line.back();
} /* END fileExtension */

void PDB2Map::setPeriodicTable()
{
    // fill in periodic table
    el2mass["H"]        = 1;
    el2mass["B"]        = 10;
    el2mass["C"]        = 12;
    el2mass["N"]        = 14;
    el2mass["O"]        = 16;
    el2mass["F"]        = 19;
    el2mass["Al"]       = 27;
    el2mass["Si"]       = 28;
    el2mass["P"]        = 31;
    el2mass["S"]        = 32;
    el2mass["Cl"]       = 35;
    el2mass["Ir"]       = 192;
    return;
} /* END setPeriodicTable */

void PDB2Map::isConvertable()
{
    string element;
    vector<Atom*> atoms = _MDtop.Atoms();
    vector<Atom*>::iterator atom_it;
    for (atom_it = atoms.begin(); atom_it != atoms.end(); atom_it++)
    {
        element = (*atom_it)->getElement();
        ba::trim(element);
        if (element.empty())
            error1("... ... PDB isn't good for conversion to QM.\n"
                   "... ... Add chemical elements or make it XYZ.");
    }
    cout << "\n... ... PDB can be converted. Continue.";
    return;
} /* END isConvertable */

} /* END namespace ctp */ } /* END namespace votca */


#endif
