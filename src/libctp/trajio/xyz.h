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
#include <boost/format.hpp>

namespace votca { namespace ctp {

//      Note:
// <xyz> is in : Angstrom
// VOTCA is in : nm
// conversion factor ( x 0.1 )
    
class xyzTrajIO: public trajIO
{
public:
    xyzTrajIO(){}
    ~xyzTrajIO(){}
    void read( const std::string&,Topology*);
    void write(const std::string&,Topology*);
private:
    void Erorr_message(string line){ // formatted and visible error message
        cout<<"\n\n"; throw runtime_error(line+'\n');
    }

};

void xyzTrajIO::read(const std::string& _filename, Topology * _top)
{
    cout << flush;
    cout << "... ... Reading structure file <xyz>" << endl;
    
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
    
               // fragment
    Fragment * _fra  = _top->AddFragment("fra_0");
                       _mol->AddFragment(_fra);
                       _seg->AddFragment(_fra);
               // inverse
               _fra->setTopology(_top);
               _fra->setMolecule(_mol);
               _fra->setSegment(_seg );
    
               
    // try to read <xyz> file
    std::ifstream _file( _filename.c_str());
    if (_file.fail())
        Erorr_message(  "... ... Can't open: " + _filename + "\n"
                        "... ... Does it exist?               \n"   
                        "... ... Is it correct file name?     \n");
    else
        cout << "... ... From: " << _filename << endl;
    
    
    // read <xyz> line by line
    string _line;
    
    // <xyz>: first two lines are tech specs
    try
    {   
        // first line, number of atoms in <xyz>
        // check if converts to int
        std::getline(_file, _line,'\n');
        boost::algorithm::trim(_line);
        int natoms = boost::lexical_cast<double>(_line);
    }
    catch(boost::bad_lexical_cast &)
    {
        Erorr_message(  "... ... Bad <xyz> file format!            \n"
                        "... ... First line : total number of atoms\n"
                        "... ... Second line: comment line         \n"
                        "... ... Broken line looks like this:      \n\n"
                        "" + _line );
    }
    
    // sample of the <xyz> line
    std::string xyz_sample = " X   -1.000   -1.000   -1.000 \n";
    
    // second line, comment line, ignore
    std::getline(_file, _line,'\n');
    
    while ( std::getline(_file, _line,'\n') ){

        // tokenize wrt space (free format)
        Tokenizer xyztok( _line, " ");
        vector<string> xyzline;
        xyztok.ToVector(xyzline);
        
        if (xyzline.size()!=4)
            Erorr_message("... ... Bad <xyz> file format!       \n"
                          "... ... Format sample:               \n\n"
                          "" + xyz_sample +                    "\n"
                          "... ... Broken line looks like this: \n\n"
                          "" + _line );
        
        string _atName     (xyzline[0]); // str,  Atom name
        string _x          (xyzline[1]); // 
        string _y          (xyzline[2]); // 
        string _z          (xyzline[3]); // 
        
        // try to convert <xyz> coords to double
        double _xd(0),_yd(0),_zd(0);
        try{
            _xd = boost::lexical_cast<double>(_x);
            _yd = boost::lexical_cast<double>(_y);
            _zd = boost::lexical_cast<double>(_z);
        }
        catch(boost::bad_lexical_cast &){
            Erorr_message("... ... Bad <xyz> file format!       \n"
                          "... ... Can't convert <xyz> string to numbers! \n"
                          "... ... Format sample:               \n\n"
                          "" + xyz_sample +                    "\n"
                          "... ... Broken line looks like this: \n\n"
                          "" + _line ); }
        
        // convert Angstroms <xyz> to nm Votca/GROMACS
        vec r(_xd , _yd , _zd); // was in Angstroms
        r = r * 0.1;            // now in nanometers
        
        // set atom
        // reconnect to topology, molecule, segment, fragment
        Atom * _atom  = _top->AddAtom(_atName);
                        _mol->AddAtom(_atom);
                        _seg->AddAtom(_atom);
                        _fra->AddAtom(_atom);
                // inverse
                _atom->setTopology(_top);
                _atom->setMolecule(_mol);        
                _atom->setSegment( _seg);
                _atom->setFragment(_fra);
        
        // set atom name, position
        _atom->setElement(_atName);
        _atom->setPos(r);
        
    } /* END of while loop */
    
    
    _file.close();
    
    return;
}

void xyzTrajIO::write(const std::string& _filename, Topology * _top)
{
    using boost::format;
    
    // stream to write on
    std::stringstream ss;
    
    // first line, atom number
    ss << boost::format("%1%\n") % _top->Atoms().size();
    
    // second line, random comment
    ss << "comment line\n";
    
    // molecule structure
    // openbabel <xyz> format: E X Y Z
    vector<Atom*>::iterator ita;
    for ( ita = _top->Atoms().begin(); ita != _top->Atoms().end(); ++ita ){
        
        format frmt("%|3|%|+10.5f|%|+10.5f|%|+10.5f|\n");
        frmt % (*ita)->getElement();
        frmt % float( (*ita)->getPos().getX() * 10. );
        frmt % float( (*ita)->getPos().getY() * 10. );
        frmt % float( (*ita)->getPos().getZ() * 10. );
        ss << frmt;
    }
    
    // write to file in <xyz> format
    // open, dump to file and close
    std::ofstream _file( _filename.c_str());
    _file << ss.rdbuf();
    _file.close();
    
    return ;    
}

} /*namespace votca END */ } /* namespace ctp END */
#endif /* VOTCA_CTP_XYZ_TRAJ_IO_H */
