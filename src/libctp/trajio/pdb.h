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
    void read(std::string,Topology*);
    void write();
    void error1(string line){ cout << endl; throw runtime_error(line); }
private:
    
};

void pdbTrajIO::read(std::string _filename, Topology* _topPtr)
{
// set molecule >> segment >> fragment
// reconnect them all
    
// molecule
    Molecule * _molPtr = 0;
    // direct
    _molPtr = _topPtr->AddMolecule("M1");
              // inverse
              _molPtr->setTopology(_topPtr);
    
// segment
    Segment  * _segPtr  = 0;
    // direct
    _segPtr = _topPtr->AddSegment("S1");
              _molPtr->AddSegment(_segPtr);
              // inverse
              _segPtr->setTopology(_topPtr);
              _segPtr->setMolecule(_molPtr);

// try to read PDB file
    std::ifstream _file( _filename.c_str());
    if (_file.fail()) {
        error1( "... ... Can not open: " + _filename + "\n"
                "... ... Does it exist? Is it correct file name?\n");
    }
    else{
        cout << endl << 
                ("... ... File opened: " + _filename );
    }

// read PDB line by line
    string _line;
    
// counters for loops
    int _newResNum = 0;
//    bool warning_showed = false;

    while ( std::getline(_file, _line,'\n') ){
        if(     boost::find_first(_line, "ATOM"  )   || 
                boost::find_first(_line, "HETATM")      
                ){
            
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
                error1( "... ... Can not convert PDB coord line!\n"
                        "... ... Error at atom index: " + _atNum + "\n"
                        "... ... Make sure this line is PDB correct\n");
            }
            
            // conversion!
            vec r(_xd , _yd , _zd); // in Angstroms
            r = r * 0.1;            // now nm

            // set fragment
            // reconnect to topology, molecule, segment
            Fragment * _fragPtr = 0;
            // make new frag for new res number
            // otherwise use last created
            if ( _newResNum != _resNumInt ){

                _newResNum = _resNumInt;
                string _newResName = _resName+'_'+_resNum;
                
                // direct
                _fragPtr = _topPtr->AddFragment(_newResName);
                           _molPtr->AddFragment(_fragPtr);
                           _segPtr->AddFragment(_fragPtr);
                          // inverse
                          _fragPtr->setTopology(_topPtr);
                          _fragPtr->setMolecule(_molPtr);
                          _fragPtr->setSegment(_segPtr);        
            }
            else{
                _fragPtr = _topPtr->Fragments().back();
            }
            if (_fragPtr==0) {error1("Zero pointer in GRO reader. Why?");}
                        
            // set atom
            // reconnect to topology, molecule, segment, fragment
            Atom * _atmPtr = 0;
            // direct
            _atmPtr = _topPtr->AddAtom(_atName);
                      _molPtr->AddAtom(_atmPtr);
                      _segPtr->AddAtom(_atmPtr);
                     _fragPtr->AddAtom(_atmPtr);
                      // inverse
                      _atmPtr->setTopology(_topPtr);
                      _atmPtr->setMolecule(_molPtr);        
                      _atmPtr->setSegment(_segPtr);
                      _atmPtr->setFragment(_fragPtr);
                      
            _atmPtr->setResnr        (_resNumInt);
            _atmPtr->setResname      (_resName);
            _atmPtr->setPos          (r);
            _atmPtr->setElement      (_atElement);
        }
    }    
    
    return;
}

void pdbTrajIO::write()
{
    std::cout << std::endl << "\n... ... I write PDB file";
}

} /*namespace votca END */ } /* namespace ctp END */
#endif /* VOTCA_CTP_PDB_TRAJ_IO_H */
