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

#ifndef VOTCA_CTP_XYZ_TRAJ_IO_H
#define VOTCA_CTP_XYZ_TRAJ_IO_H

#include <votca/ctp/trajio.h>
#include <boost/algorithm/string.hpp>

namespace votca { namespace ctp {
    namespace ba = boost::algorithm;
//      Note:
// XYZ is in   : Angstrom
// VOTCA is in : nm
// conversion factor ( x 0.1 )
    
class xyzTrajIO: public trajIO
{
public:
    xyzTrajIO(){}
    ~xyzTrajIO(){}
    void read(std::string,Topology*);
    void write();
    void error1(string line){ cout << endl; throw runtime_error(line); }
private:
    
};

void xyzTrajIO::read(std::string _filename, Topology* _topPtr)
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

// i dont create Fragment here
// because xyz has no residues
// they will be set for mapping further
                
//    Fragment * _fragPtr = 0;
//    // direct
//    _fragPtr = _topPtr->AddFragment("F1");
//               _molPtr->AddFragment(_fragPtr);
//               _segPtr->AddFragment(_fragPtr);
//                // inverse
//              _fragPtr->setTopology(_topPtr);
//              _fragPtr->setMolecule(_molPtr);
//              _fragPtr->setSegment(_segPtr);
    
// try to read xyz file
    std::ifstream _file( _filename.c_str());
    if (_file.fail()) {
        error1( "... ... Can not open: " + _filename + "\n"
                "... ... Does it exist? Is it correct file name?\n");
    }
    else{
        cout << endl << 
                ("... ... File opened: " + _filename );
    }
    
// read XYZ line by line
    string _line;
    
// XYZ: first two lines are tech specs
// XYZ check: if first line can cast to int, then ok
    try
    {   
        // first line, number of atoms in XYZ
        std::getline(_file, _line,'\n');
        ba::trim(_line);
        int numXYZatoms = boost::lexical_cast<double>(_line);
    }
    catch(boost::bad_lexical_cast &)
    {
        error1( "... ... Bad XYZ file format!\n"
                "... ... Can't find number of atoms in first line\n"
                "... ... First line: total number of atoms\n"
                "... ... Second line: comments\n"
                "... ... The line is: \n\n"
                "" + _line + "\n\n");
    }
    
    // ignore second line, it's a comment
    std::getline(_file, _line,'\n');
    
    while ( std::getline(_file, _line,'\n') ){

        // tokenize wrt space (free format)
        Tokenizer tokLine( _line, " ");
        vector<string> vecLine;
        tokLine.ToVector(vecLine);
        
        if (vecLine.size()!=4){
            error1("... ... Bad coord line in XYZ. Fix your XYZ file!\n"
                         "... ... I don't like this string: \n\n"
                         "" + _line + "\n\n"
                         "... ... Check if coords are numbers!\n");
        }
        
        string _atName     (vecLine[0]); // str,  Atom name
        string _x          (vecLine[1]); // 
        string _y          (vecLine[2]); // 
        string _z          (vecLine[3]); // 
        
        // try transform xyz coords to double
        double _xd(0),_yd(0),_zd(0);
        try{
            _xd = boost::lexical_cast<double>(_x);
            _yd = boost::lexical_cast<double>(_y);
            _zd = boost::lexical_cast<double>(_z);
        }
        catch(boost::bad_lexical_cast &)
        {
                 error1( "... ... Can't make numbers from strings.\n"
                         "... ... I don't like this string: \n\n"
                         "" + _line + "\n\n"
                         "... ... Check if coords are numbers!\n");
        }
        // conversion!
        vec r(_xd , _yd , _zd); // in Angstroms
        r = r * 0.1;            // in nanometers
        
        // set atom
        // reconnect to topology, molecule, segment, fragment
        Atom * _atmPtr = 0;
        // direct
        _atmPtr = _topPtr->AddAtom(_atName);
                _molPtr->AddAtom(_atmPtr);
                 _segPtr->AddAtom(_atmPtr);
//                _fragPtr->AddAtom(_atmPtr);
                    // inverse
                    _atmPtr->setTopology(_topPtr);
                    _atmPtr->setMolecule(_molPtr);        
                    _atmPtr->setSegment(_segPtr);
//                    _atmPtr->setFragment(_fragPtr);
        
        // set atom name, position
        _atmPtr->setElement(_atName);
        _atmPtr->setPos(r);
    }

    return;
}

void xyzTrajIO::write()
{
    std::cout << std::endl << "\n... ... I write PDB file";
}

} /*namespace votca END */ } /* namespace ctp END */
#endif /* VOTCA_CTP_XYZ_TRAJ_IO_H */
