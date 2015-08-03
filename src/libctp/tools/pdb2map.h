#ifndef _PDB2Map_H
#define _PDB2Map_H


#include <votca/ctp/topology.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/logger.h>
#include <boost/algorithm/string.hpp>
#include <votca/tools/vec.h>
#include <boost/format.hpp>
#include <algorithm>

#include <votca/ctp/trajiofactory.h>

namespace votca { namespace ctp {
    namespace ba = boost::algorithm;
//    using namespace std;

class map_vec
{
public:
    std::vector<std::string*> v;
};


    
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
    
    
private:

    string      _input_md_file;
    string      _input_qm_file;
    string      _output_file;
    
    Topology    _MDtop;
    Topology    _QMtop;

    // element:mass map
    map <string,int> el2mass;

    // internal methods
    void fileExtension(const string & name, string & ext){ 
        // return file extension
        votca::tools::Tokenizer tok(name,".");
        vector<string> factor;
        tok.ToVector(factor);
        ext = factor.back();
        return ; 
    }
    void setPeriodicTable(); // fill periodic table -> elemNumber:mass
    void Erorr_message(string line){ // formatted and visible error message
        cout<<"\n\n"; throw runtime_error(line+'\n');
    }
    
    template<class T>
    bool hasElement(vector<T>& v, T x){
        return std::find(v.begin(), v.end(), x) != v.end();
    }
    
    // formats top to <mdatoms> :               requires res
    std::string atom2map_md(Atom&);    
    // formats top to <qmatoms> and <mpoles> :  requires element
    std::string atom2map_qm(Atom&);    
    
    bool ifResChanged(Atom& a1){
        // checks if residue has chaged
        bool outcome = false;
        
        if (a1.getId() == 1)
            // first, change by default
            outcome = true;
        else
        {
          if ( a1.getFragment()->getName() !=
            a1.getTopology()->getAtom( a1.getId()-1 )->getFragment()->getName())
              // if previous is not the same, it changed
              outcome = true;
        }
        
        return outcome;
    }

};

void PDB2Map::Initialize(Property* options) 
{   
    // skip init
    cout << endl;

    // read options    
    string key = "options.pdb2map.";
    
    // md file
    if ( options->exists(key+"md") )
        _input_md_file      = options->get(key+"md").as<string> ();
    else 
        Erorr_message("... ... Option <md> not found; I need it.");
    
    // qm file
    if ( options->exists(key+"qm") )
        _input_qm_file      = options->get(key+"qm").as<string> ();
    else 
        Erorr_message("... ... Option <qm> not found; I need it.");
    
    // output xml file
    if ( options->exists(key+"output") )
        _output_file      = options->get(key+"output").as<string> ();
    else
        _output_file = "system_output.xml"; 
    
    // register all IOtraj formats
    // currently includes: pdb, gro, xyz
    trajIOFactory::RegisterAll();
    
    // show what files have been used
    cout <<    "... ... Using the following files... "   << endl <<
               "... ... MD input:\t"   << _input_md_file << endl <<
               "... ... QM input:\t"   << _input_qm_file << endl <<
               "... ... XML output:\t" << _output_file   << flush;
}

bool PDB2Map::Evaluate() {
    
    // skip 
    cout << endl;
    
    // get extensions
    string md_ext, qm_ext;
    fileExtension(_input_md_file, md_ext);
    fileExtension(_input_qm_file, qm_ext);
    
    // viable options for mapping
    vector<string> _possible_md_ext;
    _possible_md_ext.push_back("pdb");
    _possible_md_ext.push_back("gro");
    
    vector<string> _possible_qm_ext;
    _possible_qm_ext.push_back("pdb");
    _possible_qm_ext.push_back("xyz");
    
    // check if in viable options before merging into map
    if ( !hasElement<string>(_possible_md_ext, md_ext) )
        Erorr_message("... ... Input <md> type error  \n"
            "... ... Can't use: " + _input_md_file + "\n"
            "... ... Options are: pdb, gro            \n");
    if ( !hasElement<string>(_possible_qm_ext, qm_ext) )
        Erorr_message("... ... Input <qm> type error  \n"
            "... ... Can't use: " + _input_qm_file + "\n"
            "... ... Options are: pdb, xyz            \n");
    
    // read in via trajIO factory
    trajIO * 
    trajio = trajIOs().Create( md_ext );
    trajio->read(_input_md_file, &_MDtop);
    
    trajio = trajIOs().Create( qm_ext );
    trajio->read(_input_qm_file, &_QMtop);
    
    adaptQM2MD();
    
//    cout << endl << "... ... QM topology";
//    
//    fileExtension(extension,_input_qm_file);
//    trajPtr = trajIOs().Create(extension);
//    trajPtr->read(_input_qm_file,&_QMtop);
//    
//    // if "can convert", try to convert
//    // else XYZ is supplied, adapt it to MD
//    if (_conversion_possible) isConvertable();
//    else adaptQM2MD();
//    
//    topMdQm2xml();
//    
//////    LOG( logINFO, _log ) << "Reading from: " << _input_file << flush;    
//////    cout << _log;
    return true;
}

void PDB2Map::adaptQM2MD()
{
    // check if atom number is the same in QM and MD
    int numMDatoms = _MDtop.getMolecule(1)->NumberOfAtoms();
    int numQMatoms = _QMtop.getMolecule(1)->NumberOfAtoms();
    
    bool _QM2MDcompatible = (numMDatoms == numQMatoms) ? true : false;
    
    if (!_QM2MDcompatible)
        Erorr_message("... ... QM and MD molecules are not the same. \n"
                      "... ... Abort! \n");
    
//    int Natoms = _MDtop.Atoms().size();
    
//    vector<std::string> md_strings;
//    for ( itmd = _MDtop.Atoms().begin(); itmd != _MDtop.Atoms().end(); ++itmd )
//        md_strings.push_back( atom2map_md( *(*itmd) ) );
//    
//    vector<std::string> qm_strings;
//    for ( itqm = _QMtop.Atoms().begin(); itqm != _QMtop.Atoms().end(); ++itqm )
//        qm_strings.push_back( atom2map_qm( *(*itqm) ) );

    
    for ( vector<Fragment*>::iterator 
                itmdfra =  _MDtop.Fragments().begin();
                itmdfra != _MDtop.Fragments().end();
                ++itmdfra       ){
        
        cout << "New Fragment" << endl;
        
        for ( vector<Atom*>::iterator 
                itmdat =  (*itmdfra)->Atoms().begin(); 
                itmdat != (*itmdfra)->Atoms().end();
                ++itmdat        ){
        
            cout << atom2map_md( *(*itmdat) ) << '\t';
            
        }
    }
    
    
    
////        cout << atom2map_md( *(*itmd) ) << '\t' ;
//        cout << (*itmd)->getFragment()->getName();
////        if ( itmd != _MDtop.Atoms().begin() || (*itmd)->getFragment()->getName() != (*(itmd-1))->getFragment()->getName() ) 
////            cout << '\n';
//    }
    
    
    
//    // if so, proceed
//    if (_QM2MDcompatible){
//        Molecule * MDmolecule = _MDtop.getMolecule(1);
//        Molecule * QMmolecule = _QMtop.getMolecule(1);
//        
//        vector<Segment*> MDsegments = MDmolecule->Segments();
//        vector<Segment*> QMsegments = QMmolecule->Segments();
//        
//        vector<Segment*>::iterator MDSegIt;
//        vector<Segment*>::iterator QMSegIt;
//        
//        for(MDSegIt = MDsegments.begin(), QMSegIt = QMsegments.begin();
//                MDSegIt < MDsegments.end();
//                MDSegIt++, QMSegIt++  )
//        {
//            Fragment * QMfragment = 0;
//            
//            vector<Atom*> MDSegAtoms = (*MDSegIt)->Atoms();
//            vector<Atom*> QMSegAtoms = (*QMSegIt)->Atoms();
//            
//            vector<Atom*>::iterator MDSegAtIt;
//            vector<Atom*>::iterator QMSegAtIt;
//            
//            int old_res_num(-1),new_res_num(-1);
//            string res_name = "bad_wolf";
//            
//            for(MDSegAtIt = MDSegAtoms.begin(), QMSegAtIt = QMSegAtoms.begin();
//                    MDSegAtIt < MDSegAtoms.end();
//                    MDSegAtIt++, QMSegAtIt++)
//            {
//                new_res_num = (*MDSegAtIt)->getResnr();
//                if ( new_res_num != old_res_num)
//                {
//                    old_res_num = new_res_num;
//                    res_name = (*MDSegAtIt)->getResname();
//                    
//                    QMfragment = _QMtop.AddFragment(res_name);
//                    QMmolecule->AddFragment(QMfragment);
//                    (*QMSegIt)->AddFragment(QMfragment);
//                    
//                    QMfragment->setTopology(&_QMtop);
//                    QMfragment->setMolecule(QMmolecule);
//                    QMfragment->setSegment(*QMSegIt);
//                    
//                    (*QMSegAtIt)->setFragment(QMfragment);
//                    QMfragment->AddAtom(*QMSegAtIt);                    
//                }
//                else
//                {
//                    (*QMSegAtIt)->setFragment(QMfragment);
//                    QMfragment->AddAtom(*QMSegAtIt);    
//                }
//            }
//        }
//    }
//    else{
//        error1("... ... Number of MD atoms is different from QM.\n"
//               "... ... If it's the case of reduced molecule, "
//                        " I need a map.\n"
//               " ... ... Tags: map\n"
//               " ... ... NOT IMPLEMENTED\n");
//    }
} /* END adaptQM2MD */

void PDB2Map::topMdQm2xml(){return;}
//            if ( (!_has_xyz && !warning_showed) || _line.size()<78 ){
//            if ( (!_has_xyz && !warning_showed) || _line.size()<78 ){
//            if ( (!_has_xyz && !warning_showed) || _line.size()<78 ){
//            if ( (!_has_xyz && !warning_showed) || _line.size()<78 ){
//               cout << endl << "... ... WARNING: No chemical elements in PDB!\n"
//                            << "... ... Expect: empty slots \n"
//                            << "in <qmatoms> and <multipoles>, "
//                               "zeros in <weights>.\n"
//                            << "... ... To add chemical symbols use: "
//                               "editconf (GROMACS), babel, "
//                               "(hands+pdb format)";                   
//               warning_showed = true;
//            }
//            std::cout << "The size of str is " << _line.size() << " bytes.\n";
//            std::cout << "The size of str is " << _line.size() << " bytes.\n";
//            std::cout << "The size of str is " << _line.size() << " bytes.\n";
//            std::cout << "The size of str is " << _line.size() << " bytes.\n";
//{
//    cout << endl << "... ... (A)merging XML from MD and QM topologies.";
//
//    Molecule * MDmolecule = _MDtop.getMolecule(1);
//    Molecule * QMmolecule = _QMtop.getMolecule(1);
//    
//    //
//    // xml stuff
//    //
//    
//    Property record;
//    Property *ptopology_p = &record.add("topology","");
//    Property *pmolecules_p = &ptopology_p->add("molecules","");
//    Property *pmolecule_p = &pmolecules_p->add("molecule","");
//    pmolecule_p->add("name","random_molecule_name");
//    pmolecule_p->add("mdname","name_from_topol.top");
//    Property *psegments_p = &pmolecule_p->add("segments","");
//    Property *psegment_p = &psegments_p->add("segment","");
//    psegment_p->add("name","random_segment_name");
//    
//    // qc data
//    psegment_p->add("qmcoords","QC_FILES/your_file_with.xyz");
//    psegment_p->add("orbitals","QC_FILES/your_file_with.fort7");
//    psegment_p->add("basisset","INDO");
//    
//    // hole
//    psegment_p->add("torbital_h","number_of_electrons,150");
//    psegment_p->add("U_cC_nN_h","0.0000");
//    psegment_p->add("U_nC_nN_h","0.0000");
//    psegment_p->add("U_cN_cC_h","0.0000");
//    
//    // electron
//    psegment_p->add("torbital_e","number_of_electrons+1,151");
//    psegment_p->add("U_cC_nN_e","0.0000");
//    psegment_p->add("U_nC_nN_e","0.0000");
//    psegment_p->add("U_cN_cC_e","0.0000");
//
//    // mps entry
//    psegment_p->add("multipoles_n","MP_FILES/file_with.mps");
//    psegment_p->add("multipoles_h","MP_FILES/file_with.mps");
//    psegment_p->add("multipoles_e","MP_FILES/file_with.mps");
//    
//    psegment_p->add("map2md","0");
//    
//    // main body 
//    Property *pfragments_p = &psegment_p->add("fragments","");
//                                        
//    vector < Segment * > allMdSegments = MDmolecule->Segments();
//    vector < Segment * > allQmSegments = QMmolecule->Segments();
//  
//    vector < Segment * >::iterator segMdIt;
//    vector < Segment * >::iterator segQmIt;
//                                        
//    for ( segMdIt = allMdSegments.begin(), 
//                segQmIt = allQmSegments.begin();
//            
//          (allMdSegments.size() > allQmSegments.size()) ? 
//              segMdIt < allMdSegments.end() :
//              segQmIt < allQmSegments.end();
//            
//          segMdIt++, segQmIt++)
//    {
//        
//        vector < Fragment * > allMdFragments = (*segMdIt)->Fragments();
//        vector < Fragment * > allQmFragments = (*segQmIt)->Fragments();
//                                        
//        vector < Fragment * >::iterator fragMdIt;
//        vector < Fragment * >::iterator fragQmIt;
//                                        
//        for ( fragMdIt = allMdFragments.begin() , 
//                fragQmIt = allQmFragments.begin();
//                
//              (allMdFragments.size() > allQmFragments.size()) ?
//                  fragMdIt < allMdFragments.end() :
//                  fragQmIt < allQmFragments.end();
//                
//              fragMdIt++,fragQmIt++ )
//        {
//            string mapName;            
//            stringstream mapMdAtoms;
//            stringstream mapQmAtoms;
//            stringstream mapMpoles;
//            stringstream mapWeight;
//            stringstream mapFrame;
//
//            mapName      = (*fragMdIt)->getName() ;
//            
//            int localCounter = 0;
//            vector < Atom * > allMdAtoms = (*fragMdIt)->Atoms();
//            vector < Atom * > allQmAtoms = (*fragQmIt)->Atoms();
//
//            vector < Atom * >::iterator atomMdIt;
//            vector < Atom * >::iterator atomQmIt;
//            for ( atomMdIt = allMdAtoms.begin(),
//                    atomQmIt = allQmAtoms.begin();
//                    
//                  (allMdAtoms.size() > allQmAtoms.size()) ? 
//                      atomMdIt < allMdAtoms.end() :
//                      atomQmIt < allQmAtoms.end();
//                  
//                  atomMdIt++, atomQmIt++ )
//            {
//                if (atomMdIt < allMdAtoms.end())
//                {
//                        mapMdAtoms << boost::format("%=13s") % 
//                                (boost::format("%s:%s:%s") 
//                                   % (*atomMdIt)->getResnr()
//                                   % (*atomMdIt)->getResname()
//                                   % (*atomMdIt)->getName()).str() ;
//                        
//                                if (el2mass.find((*atomQmIt)->getElement()) 
//                                        != el2mass.end())
//                                {
//                        mapWeight << boost::format("%=13i")
//                                           % el2mass[(*atomQmIt)->getElement()];
//                                }
//                                else
//                                {
//                        mapWeight << boost::format("%=13i") % " ";        
//                                }
//                }
//                else
//                {
//                        mapMdAtoms << boost::format("%=13s") %
//                                (boost::format("%s:%s:%s") 
//                                        % " " % " " % " " ).str();
//                        
//                        mapWeight << boost::format("%=13i") % " " ;
//                }
//                
//                if (atomQmIt < allQmAtoms.end())
//                {
//                        mapQmAtoms << boost::format("%=13s") %
//                                (boost::format("%1%:%2% ") 
//                                   % (*atomQmIt)->getId()
//                                   % (*atomQmIt)->getElement()).str() ;
//                }
//                else
//                {
//                        mapQmAtoms << boost::format("%=13s") %
//                                (boost::format("%1%:%2% ") 
//                                        % " " % " " ).str();
//                }
//                
//                if (localCounter < 3 && localCounter < allMdAtoms.size())
//                {
//                        mapFrame << boost::format("%=5i")
//                                   % (*atomMdIt)->getId();
//                }
//                localCounter++;
//            }
//            
//            mapMpoles << " " << mapQmAtoms.str();
//            
//            Property *pfragment_p  = &pfragments_p->add("fragment","");
//            pfragment_p->add("name", mapName);
//            pfragment_p->add("mdatoms", mapMdAtoms.str());
//            pfragment_p->add("qmatoms", mapQmAtoms.str());
//            pfragment_p->add("mpoles",  mapMpoles.str());
//            pfragment_p->add("weights", mapWeight.str());
//            pfragment_p->add("localframe", mapFrame.str());
//
//         }
//    }
//    
////    cout << endl << setlevel(1) << XML << record;
//        
//    ofstream outfile( _output_file.c_str() );
//    
//    PropertyIOManipulator XML(PropertyIOManipulator::XML,1,""); 
//    outfile << XML << record;
//    outfile.close();
//    
//    cout << endl << "... ... Topology written to: " + _output_file;
//    cout << endl << "... ... PDB2MAP is done. Victory.";
//    
//    return;
//} /* END topMdQm2xml */

//void PDB2Map::fileExtension(string & extension, const string & filename)
//{
//    Tokenizer tokenized_line( filename, ".");
//    vector<string> vector_line;
//    tokenized_line.ToVector(vector_line);
//    extension = vector_line.back();
//} /* END fileExtension */

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

    el2mass["Ga"]       = 70;
    el2mass["Ge"]       = 73;
    el2mass["As"]       = 75;
    el2mass["Se"]       = 79;
    el2mass["Br"]       = 80;

    el2mass["Ir"]       = 192;

    return;
} /* END setPeriodicTable */

//void PDB2Map::isConvertable()
//{
//    string element;
//    vector<Atom*> atoms = _MDtop.Atoms();
//    vector<Atom*>::iterator atom_it;
//    for (atom_it = atoms.begin(); atom_it != atoms.end(); atom_it++)
//    {
//        element = (*atom_it)->getElement();
//        ba::trim(element);
//        if (element.empty())
//            error1("... ... PDB isn't good for conversion to QM.\n"
//                   "... ... Add chemical elements or make it XYZ.");
//    }
//    cout << "\n... ... PDB can be converted. Continue.";
//    return;
//} /* END isConvertable */

std::string PDB2Map::atom2map_md(Atom& atom){
    // builds blocks for <mdatoms>
    std::stringstream ss;
    std::string res_id, res_name, atom_name;
    
    Tokenizer tok( atom.getFragment()->getName(), "_"); // res: name_id
    vector<std::string> v;
    tok.ToVector(v);
    
    res_id    = v.back();  // res id
    res_name  = v.front(); // res name
    
    ss << res_id << ":" << res_name << ":" << atom.getName();
    return ss.str();
}

std::string PDB2Map::atom2map_qm(Atom& atom){
    // builds blocks for <qmatom>
    std::stringstream ss;
    ss << atom.getId() << ":" << atom.getElement();
    return ss.str();
}

} /* END namespace ctp */ } /* END namespace votca */


#endif
