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
#include <boost/format.hpp>

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
    void read( const std::string&,Topology*);
    void write(const std::string&,Topology*);
private:
    void Erorr_message(string line){ // formatted and visible error message
        cout<<"\n\n"; throw runtime_error(line+'\n');
    }    
};

void groTrajIO::read(const std::string& _filename, Topology* _top)
{
    cout << flush;
    cout << "... ... Reading structure file <gro>" << endl;
    
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
    

    // try to read <gro> file
    std::ifstream _file( _filename.c_str());
    if (_file.fail())
        Erorr_message(  "... ... Can't open: " + _filename + "\n"
                        "... ... Does it exist?               \n"   
                        "... ... Is it correct file name?     \n");
    else
        cout << "... ... From: " << _filename << endl;
    
    
    // read <gro> line by line
    string _line;
    int _natoms = 0;  // total num of atoms in GRO
    
    // <gro>: first two lines are tech specs
    // first line, comment line, skip
    std::getline(_file, _line,'\n');

    // <gro>: second line, number of atoms
    // check if convert to int
    try
    {   
        std::getline(_file, _line,'\n');
        ba::trim(_line);
        _natoms = boost::lexical_cast<int>(_line);
    }
    catch(boost::bad_lexical_cast &)
    {
        Erorr_message(  "... ... Bad <gro> file format!            \n"
                        "... ... First line : comment line         \n"
                        "... ... Second line: total number of atoms\n"
                        "... ... Broken line looks like this:      \n\n"
                        "" + _line );
    }
    
    // sample of the <gro> line
    std::string gro_sample = "    1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434 \n";
    
    // loop over <gro> atoms
    for (int _atCounter=0; _atCounter < _natoms; ++_atCounter ){
            
        std::getline(_file, _line,'\n');
            
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
          Erorr_message("... ... Bad <gro> file format!       \n"
                        "... ... Can't convert <gro> string to numbers! \n"
                        "... ... Format sample:               \n\n"
                        "" + gro_sample +                    "\n"
                        "... ... Broken line looks like this: \n\n"
                        "" + _line );
        }

        vec r(_xd , _yd , _zd);
        
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
        
    } /* END of <gro> atoms loop */
    
    // last line, box dimensions
    std::getline(_file, _line,'\n');
    
    // is this a correct format ?
    _top->setBox(votca::tools::matrix(vec(_line),vec(0,0,0),vec(0,0,0)));

    
    _file.close();
    
    return;
}

void groTrajIO::write(const std::string& _filename, Topology * _top)
{
    using boost::format;
    
    // stream to write on
    std::stringstream ss;
    
    // first line, comment
    ss << "GRO file::Generated by VOTCA-CTP \n";
    
    // second line, atom number
    ss << boost::format("%1%\n") % _top->Atoms().size();
    
    // molecule structure
    // GROMACS <GRO> format
    vector<Atom*>::iterator ita;
    for ( ita = _top->Atoms().begin(); ita != _top->Atoms().end(); ++ita ){
        
        int         res_id   = (*ita)->getFragment()->getId();
        std::string res_name = *(Tokenizer((*ita)->getFragment()->getName(), "_").begin());
        std::string at_name  = (*ita)->getName();
        int         at_id    = (*ita)->getId();
        
        format frmt("%|5|%|5|%|5|%|5|%|+8.3f|%|+8.3f|%|+8.3f|%|+8.4f|%|+8.4f|%|+8.4f|\n");
        frmt % res_id;
        frmt % res_name;
        frmt % at_name;
        frmt % at_id;
        // positions
        frmt % (*ita)->getPos().getX();
        frmt % (*ita)->getPos().getY();
        frmt % (*ita)->getPos().getZ();
        // velocities
        frmt % 0.0;
        frmt % 0.0;
        frmt % 0.0;

        ss << frmt;
    }
    
    // final line, box dims
    format frmt_bx("%|+10.5f|%|+10.5f|%|+10.5f|\n");
    frmt_bx % _top->getBox().getCol(0).getX();
    frmt_bx % _top->getBox().getCol(0).getY();
    frmt_bx % _top->getBox().getCol(0).getZ();
    ss << frmt_bx;
    
    // write to file in <gro> format
    // open, dump to file and close
    std::ofstream _file( _filename.c_str());
    _file << ss.rdbuf();
    _file.close();
    
    return ;    
}

} /* namespace votca END */ } /* namespace ctp END */
#endif /* VOTCA_CTP_GRO_TRAJ_IO_H */
