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

#ifndef VOTCA_CTP_PDB_TRAJ_IO_H
#define VOTCA_CTP_PDB_TRAJ_IO_H

#include <votca/ctp/trajio.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

namespace votca { namespace ctp {
    namespace ba = boost::algorithm;
//      Note:
// PDB is in   : Angstrom
// VOTCA is in : nm
// conversion factor ( x 0.1 )    
    
class pdbTrajIO: public trajIO
{
public:
    pdbTrajIO(){}
    ~pdbTrajIO(){}
    void read( const std::string&,Topology*);
    void write(const std::string&,Topology*);
private:
    void Erorr_message(string line){ // formatted and visible error message
        cout<<"\n\n"; throw runtime_error(line+'\n');
    }
};

void pdbTrajIO::read(const std::string& _filename, Topology* _top)
{
    cout << flush;
    cout << "... ... Reading structure file <pdb>" << endl;
    
    // set molecule >> segment >> fragment
    // reconnect them all

    // molecule    
    Molecule * _mol = _top->AddMolecule("mol_0");
               // inverse
               _mol->setTopology(_top);
    
    // segment    
    Segment  * _seg  = _top->AddSegment("seg_0");
                       _mol->AddSegment(_seg);
               // inverse
               _seg->setTopology(_top);
               _seg->setMolecule(_mol);
    
    // fragment, to be set up later
    Fragment * _fra = 0;
    

    // try to read <pdb> file
    std::ifstream _file( _filename.c_str());
    if (_file.fail())
        Erorr_message(  "... ... Can't open: " + _filename + "\n"
                        "... ... Does it exist?               \n"   
                        "... ... Is it correct file name?     \n");
    else
        cout << "... ... From: " << _filename << endl;
    
    std::string pdb_sample;
    pdb_sample = "ATOM     34  CA AARG A  -3      12.353  85.696  94.456  0.50 36.67           C  \n";
    
    // read <pdb> lines
    std::string _line;
    while ( std::getline(_file, _line,'\n') )
        if(     boost::find_first(_line, "ATOM"  )   || 
                boost::find_first(_line, "HETATM")      
                ){
            
            if ( _line.size() < 80 )
            {
                // if PDB string is truncated, restore
                ba::trim(_line);
                _line += std::string( 80-_line.size(), ' ');
            }
            
            //      according to PDB format
            string _recType    (_line,( 1-1),6); // str,  "ATOM", "HETATM"
            string _atNum      (_line,( 7-1),6); // int,  Atom serial number
            string _atName     (_line,(13-1),4); // str,  Atom name
            string _atAltLoc   (_line,(17-1),1); // char, Alternate location indicator
            string _resName    (_line,(18-1),4); // str,  Residue name
            string _chainID    (_line,(22-1),1); // char, Chain identifier
            string _resNum     (_line,(23-1),4); // int,  Residue sequence number
            string _atICode    (_line,(27-1),1); // char, Code for insertion of res
            string _x          (_line,(31-1),8); // float 8.3 ,x
            string _y          (_line,(39-1),8); // float 8.3 ,y
            string _z          (_line,(47-1),8); // float 8.3 ,z
            string _atOccup    (_line,(55-1),6); // float  6.2, Occupancy
            string _atTFactor  (_line,(61-1),6); // float  6.2, Temperature factor
            string _segID      (_line,(73-1),4); // str, Segment identifier
            string _atElement  (_line,(77-1),2); // str, Element symbol
            string _atCharge   (_line,(79-1),2); // str, Charge on the atom

            ba::trim(_recType);
            ba::trim(_atNum);
            ba::trim(_atName);
            ba::trim(_atAltLoc);
            ba::trim(_resName);
            ba::trim(_chainID);
            ba::trim(_resNum);
            ba::trim(_atICode);
            ba::trim(_x);
            ba::trim(_y);
            ba::trim(_z);
            ba::trim(_atOccup);
            ba::trim(_atTFactor);
            ba::trim(_segID);
            ba::trim(_atElement);
            ba::trim(_atCharge);
            
            double _xd(0),_yd(0),_zd(0);
            int _resNumInt(0); 
            
            try
            {
                _xd = boost::lexical_cast<double>(_x);
                _yd = boost::lexical_cast<double>(_y);
                _zd = boost::lexical_cast<double>(_z);
                _resNumInt = boost::lexical_cast<int>(_resNum);
            }
            catch(boost::bad_lexical_cast &)
            {
                Erorr_message("... ... Bad <pdb> file format!       \n"
                              "... ... Can't convert <gro> string to numbers! \n"
                              "... ... Format sample:               \n\n"
                              "" + pdb_sample +                    "\n"
                              "... ... Broken line looks like this: \n\n"
                              "" + _line );
            }
            
            vec r(_xd , _yd , _zd);
            r = r * 0.1 ; // angstrom to nm 
            
            string _newFragName = _resName+'_'+_resNum;

            if (_top->Fragments().empty()){
                // if empty, create first fragment
                _fra =  _top->AddFragment(_newFragName);
                        _mol->AddFragment(_fra);
                        _seg->AddFragment(_fra);
                          // inverse
                          _fra->setTopology(_top);
                          _fra->setMolecule(_mol);
                          _fra->setSegment(_seg);   
            }
            else if (_top->Fragments().back()->getName() != _newFragName){
                // if name does not match, create new fragment
                _fra =  _top->AddFragment(_newFragName);
                        _mol->AddFragment(_fra);
                        _seg->AddFragment(_fra);
                          // inverse
                          _fra->setTopology(_top);
                          _fra->setMolecule(_mol);
                          _fra->setSegment(_seg);   
            }
            
            Atom * _atom  = _top->AddAtom(_atName);
                            _mol->AddAtom(_atom);
                            _seg->AddAtom(_atom);
                            _fra->AddAtom(_atom);
                      // inverse
                      _atom->setTopology(_top);
                      _atom->setMolecule(_mol);        
                      _atom->setSegment( _seg);
                      _atom->setFragment(_fra);

            _atom->setResnr        (_resNumInt);
            _atom->setResname      (_resName);
            _atom->setPos          (r);
            _atom->setElement      (_atElement);
            
    } /*  END of while loop */    
    
    
    _file.close();
    
    return;
}

void pdbTrajIO::write(const std::string& _filename, Topology * _top)
{
    using boost::format;
    
    // stream to write on
    std::stringstream ss;
    
    // first line, comment
    ss << "TITLE     PDB file::Generated by VOTCA-CTP \n";
    
    // second line, MODEL modifier
    ss << boost::format("%|-6|    %|4|\n") % "MODEL" % 1;
    
    // molecule structure
    // GROMACS <PDB> format
    vector<Atom*>::iterator ita;
    for ( ita = _top->Atoms().begin(); ita != _top->Atoms().end(); ++ita ){
        
        int         res_id   = (*ita)->getFragment()->getId();
        std::string res_name = *(Tokenizer((*ita)->getFragment()->getName(), "_").begin());
        std::string at_name  = (*ita)->getName();
        int         at_id    = (*ita)->getId();
        
        format frmt("%|-6|%|5| %|-4|%|1|%|3| %|1|%|4|%|1|   %|8.3f|%|8.3f|%|8.3f|%|6.2f|%|6.2f|          %|2|%|2|\n");
        // marker
        frmt % "ATOM";
        // atom data
        frmt % at_id;
        frmt % at_name;
        frmt % " "; // alternate location indicator 	
        // res data 
        frmt % res_name;
        frmt % " "; // chain identifier
        frmt % res_id;
        frmt % " "; // code for insertion of residues
        // positions
        frmt % float( (*ita)->getPos().getX() * 10. );
        frmt % float( (*ita)->getPos().getY() * 10. );
        frmt % float( (*ita)->getPos().getZ() * 10. );
        // occupancy + temperature factor 
        frmt % 1.0; // occupancy
        frmt % 0.0; // temperature factor
        frmt % (*ita)->getElement(); // element symbol
        frmt % " "; // charge on the atom

        ss << frmt;
    }
    
    // final line, box dims
    ss << "TER\nENDMDL\n";
    
    // write to file in <gro> format
    // open, dump to file and close
    std::ofstream _file( _filename.c_str());
    _file << ss.rdbuf();
    _file.close();
    
    return ;    
}

} /*namespace votca END */ } /* namespace ctp END */
#endif /* VOTCA_CTP_PDB_TRAJ_IO_H */
