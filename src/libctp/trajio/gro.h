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

#ifndef VOTCA_CTP_GRO_TRAJ_IO_H
#define VOTCA_CTP_GRO_TRAJ_IO_H

#include <votca/ctp/trajio.h>
#include <boost/algorithm/string.hpp>

namespace votca { namespace ctp {
    namespace ba = boost::algorithm;
//      Note:
// PDB is in   : Angstrom
// VOTCA is in : nm
// conversion factor ( x 0.1 )    
    
class groTrajIO: public trajIO
{
public:
    groTrajIO(){}
    ~groTrajIO(){}
    void read(std::string,Topology*);
    void write();
    void error1(string line){ cout << endl; throw runtime_error(line); }
private:
    
};

void groTrajIO::read(std::string _filename, Topology* _topPtr)
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

// try to read GRO file
    std::ifstream _file( _filename.c_str());
    if (_file.fail())
        error1( "... ... Can not open: " + _filename + "\n"
                "... ... Does it exist? Is it correct file name?\n");
    else
        cout << endl << ("... ... File opened: " + _filename + "\n");

// read GRO line by line
    string _line;
    
// counters for loops
    int _newResNum = -1; // res reference
    int _atTotl = 0;  // total num of atoms in GRO
    int _atCntr = 0;  // atom number counter
    
// GRO: first two lines are tech specs -> ignore them
// ignore first line, it's a comment
    std::getline(_file, _line,'\n');

// GRO check: if second line can cast to int, then ok
    try
    {   
        // first line, number of atoms in XYZ
        std::getline(_file, _line,'\n');
        ba::trim(_line);
        _atTotl = boost::lexical_cast<int>(_line);
    }
    catch(boost::bad_lexical_cast &)
    {
        error1( "... ... Bad GRO file format!\n"
                "... ... First line contains Comment.\n"
                "... ... Second line contains Number of atoms.\n"
                "... ... I cant convert second line, check it!.\n\n"
                "" + _line + "\n\n");
    }

    // actual loop
    while ( std::getline(_file, _line,'\n') ){
        if (_atCntr < _atTotl){
            
            string _resNum     (_line, 0,5); // int,  Residue number
            string _resName    (_line, 5,5); // str,  Residue name
            string _atName     (_line,10,5); // str,  Atom name
            string _atNum      (_line,15,5); // int,  Atom number
            string _x          (_line,20,8); // float 8.3 ,x
            string _y          (_line,28,8); // float 8.3 ,y
            string _z          (_line,36,8); // float 8.3 ,z
            
            ba::trim(_atNum);
            ba::trim(_atName);
            ba::trim(_resNum);
            ba::trim(_resName);
            ba::trim(_x);
            ba::trim(_y);
            ba::trim(_z);
            
            // try cast
            int _resNumInt(0),_atNumInt(0);
            double _xd(0),_yd(0),_zd(0);
            try
            {
                _resNumInt = boost::lexical_cast<int>(_resNum);
                _atNumInt  = boost::lexical_cast<int>(_atNum);

                _xd = boost::lexical_cast<double>(_x);
                _yd = boost::lexical_cast<double>(_y);
                _zd = boost::lexical_cast<double>(_z);
            }
            catch (boost::bad_lexical_cast &)
            {
                error1( "... ... Can not convert GRO coord line!\n"
                        "... ... Atom number: " + _atNum + "\n"
                        "... ... Make sure this line is GRO style\n");
            }
            
            vec r(_xd , _yd , _zd);
                
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
        
        }
        _atCntr++;
    }

    
    
    return;
}

void groTrajIO::write()
{
    std::cout << std::endl << "\n... ... I write PDB file";
}

} /* namespace votca END */ } /* namespace ctp END */
#endif /* VOTCA_CTP_GRO_TRAJ_IO_H */
