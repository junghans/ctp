/* 
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICEN_olE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>

#include <votca/ctp/threecenters.h>

#include <votca/ctp/aobasis.h>
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>
#include <votca/ctp/logger.h>
#include <votca/tools/linalg.h>

using namespace std;
using namespace votca::tools;

namespace votca {  namespace ctp {
namespace ub = boost::numeric::ublas;

        
        void TCMatrix::Symmetrize( const ub::matrix<double>& _coulomb ){
            // cout << "lala" << endl;
            
            for ( int _i_occ = 0; _i_occ < this->get_mtot(); _i_occ++ ){
                //POTENTIALLY A BUG
                _matrix[ _i_occ ] = ub::prod(_coulomb, _matrix[ _i_occ ] );

            }
            
           
        }

        
        
        
        
        
        void TCMatrix::Fill(AOBasis& _gwbasis, AOBasis& _dftbasis, ub::matrix<double>& _dft_orbitals) {

            
            std::vector< ub::matrix<double> > _block( this->get_mtot() );
            
            
            // loop over all shells in the GW basis and get _Mmn for that shell
            for (vector< AOShell* >::iterator _is = _gwbasis.firstShell(); _is != _gwbasis.lastShell(); _is++) {
                AOShell* _shell = _gwbasis.getShell(_is);
                int _start = _shell->getStartIndex();
                //int _end = _start + _shell->getNumFunc();

       
                // each element is a shell_size-by-n matrix, initialize to zero
                for ( int i = 0; i < this->get_mtot() ; i++){
                    _block[i] = ub::zero_matrix<double>( _shell->getNumFunc() ,this->get_ntot());
                }
                
                                
                // Fill block for this shell (3-center overlap with _dft_basis + multiplication with _dft_orbitals )
                FillBlock(_block, _shell, _dftbasis, _dft_orbitals);
                
                // put into correct position

                for ( int _m_band = 0; _m_band < this->get_mtot(); _m_band++ ){
                    for ( int _i_gw = 0; _i_gw < _shell->getNumFunc(); _i_gw++ ){
                        for ( int _n_band = 0; _n_band < this->get_ntot(); _n_band++){
                            
                
                            _matrix[_m_band]( _start + _i_gw , _n_band ) = _block[_m_band]( _i_gw , _n_band );
                            
                            
                            
                        }
                    }
                }
                
                
                
            }
        }
        
        // new storage test!
        void TCMatrix::FillBlock(std::vector<  ub::matrix<double> >& _block, AOShell* _shell, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals) {

           // cout << TimeStamp() << "  ... get the matrices for GW basis block of type " << _shell->getType() << endl;

            /* calculate 3-center overlap between alpha(DFT),beta(GW), gamma(DFT)
             * - for fixed beta(GW) given in call to this function, and all 
             *    alpha,gamma combinations in DFT basis 
             * - output should be in normalized spherical Gaussians
             */

            // prepare local container
            ub::matrix<double> _imstore = ub::zero_matrix<double>(this->mtotal * _shell->getNumFunc(), dftbasis._AOBasisSize);
            
            // loop alpha(dft)
            for (vector< AOShell* >::iterator _row = dftbasis.firstShell(); _row != dftbasis.lastShell(); _row++) {
                AOShell* _shell_row = dftbasis.getShell(_row);
                int _row_start = _shell_row->getStartIndex();
                int _row_end = _row_start + _shell_row->getNumFunc();

                // get slice of _dft_orbitals for m-summation, belonging to this shell
                ub::matrix_range< ub::matrix<double> > _m_orbitals = ub::subrange(_dft_orbitals, this->mmin, this->mmax + 1, _row_start, _row_end);


                // loop gamma(dft)
                for (vector< AOShell* >::iterator _col = dftbasis.firstShell(); _col != dftbasis.lastShell(); _col++) {
                    AOShell* _shell_col = dftbasis.getShell(_col);
                    int _col_start = _shell_col->getStartIndex();
                    int _col_end = _col_start + _shell_col->getNumFunc();

                    // get 3-center overlap directly as "supermatrix"
                    ub::matrix<double> _subvector = ub::zero_matrix<double>(_shell_row->getNumFunc(), _shell->getNumFunc() * _shell_col->getNumFunc());
                    bool nonzero = FillThreeCenterOLBlock(_subvector, _shell, _shell_row, _shell_col);
          
                    // if this contributes, multiply _subvector with _dft_orbitals and place in _imstore
                    if (nonzero) {

                        ub::matrix<double> _temp = ub::prod(_m_orbitals, _subvector);
                         
                        // put _temp into _imstore
                         for ( unsigned int _i_band = 0; _i_band < _temp.size1(); _i_band++ ){
                            for ( int _i_gw = 0; _i_gw < _shell->getNumFunc() ; _i_gw++ ){
                                int _ridx = _shell->getNumFunc() * ( _i_band - this->mmin ) + _i_gw;
                                
                                for ( int _cidx = _col_start; _cidx < _col_end; _cidx++ ){
                                    
                                int _tidx = _shell_col->getNumFunc() * ( _i_gw  ) +  _cidx-_col_start ;

                                    _imstore( _ridx, _cidx ) += _temp( _i_band, _tidx ); 
                                }
                            }
                        } 
                    } // adding contribution                        
                } // gamma-loop
            } // alpha-loop


            // get transposed slice of _dft_orbitals
            ub::matrix<double> _n_orbitals = ub::trans(ub::subrange(_dft_orbitals, this->nmin, this->nmax + 1, 0, _dft_orbitals.size2()));
            
            // Now, finally multiply _imstore with _n_orbitals
            ub::matrix<double> _temp = ub::prod(_imstore, _n_orbitals);

            // and put it into the block it belongs to
            for (int _m_band = 0; _m_band < this->mtotal; _m_band++) {
                for (int _i_gw = 0; _i_gw < _shell->getNumFunc(); _i_gw++) {
                    int _midx = _shell->getNumFunc() *(_m_band - this->mmin) + _i_gw;
                    for (int _n_band = 0; _n_band < this->ntotal; _n_band++) {
                        _block[_m_band](_i_gw, _n_band) = _temp(_midx, _n_band);
//                          cout << "MNb  [" << _i_gw << ":" << _m_band << ":" << _n_band << "] =" << _block(_i_gw)( _m_band, _n_band ) << endl;
                    }
                }
            }
        }
        
        
        
        
        

 

        
        bool TCMatrix::FillThreeCenterOLBlock(ub::matrix<double>& _subvector, AOShell* _shell_gw, AOShell* _shell_alpha, AOShell* _shell_gamma) {

            // this now is supposed to determine one alpha(DFT)-beta(GW)-gamma(DFT) subblock

            // cout << "         ... trying to get the 3-center overlap block between shells of type [" << _shell_alpha->getType() << ":" << _shell_gw->getType() << ":" << _shell_gamma->getType() << "]" << endl;


            // get decay constants (this all is still valid only for uncontracted functions)
            const double& _decay_gw = (*_shell_gw->firstGaussian())->decay;
            const double& _decay_alpha = (*_shell_alpha->firstGaussian())->decay;
            const double& _decay_gamma = (*_shell_gamma->firstGaussian())->decay;

            // get shell positions
            const vec& _pos_gw = _shell_gw->getPos();
            const vec& _pos_alpha = _shell_alpha->getPos();
            const vec& _pos_gamma = _shell_gamma->getPos();

            // shell info, only lmax tells how far to go
            int _lmax_gw = _shell_gw->getLmax();
            int _lmax_alpha = _shell_alpha->getLmax();
            int _lmax_gamma = _shell_gamma->getLmax();

            // set size of internal block for recursion
            int _ngw = this->getBlockSize(_lmax_gw);
            int _nalpha = this->getBlockSize(_lmax_alpha);
            int _ngamma = this->getBlockSize(_lmax_gamma);

            // definition of cutoff for contribution
            const double gwaccuracy = 1.e-9; // should become an OPTION
            double threshold = -(_decay_alpha + _decay_gamma + _decay_gw) * log(gwaccuracy);

            // check first threshold
            vec _diff = _pos_alpha - _pos_gw;
            double test = _decay_alpha * _decay_gw * _diff*_diff;
            if (test > threshold) return false;

            // check second threshold
            _diff = _pos_gamma - _pos_gw;
            test += _decay_gamma * _decay_gw * _diff*_diff;
            if (test > threshold) return false;

            // check third threshold
            _diff = _pos_alpha - _pos_gamma;
            test += _decay_alpha * _decay_gamma * _diff*_diff;
            if (test > threshold) return false;

            // if all threshold test are passed, start evaluating

            // some helpers
            double fak = 0.5 / (_decay_alpha + _decay_gw + _decay_gamma);
            double fak2 = 2.0 * fak;
            //double fak3 = 3.0 * fak;

            vec gvv = fak2 * (_decay_alpha * _pos_alpha + _decay_gw * _pos_gw + _decay_gamma * _pos_gamma);
            vec gma = gvv - _pos_alpha;
            vec gmb = gvv - _pos_gamma;
            vec gmc = gvv - _pos_gw;

            double gma0 = gma.getX();
            double gmb0 = gmb.getX();
            double gmc0 = gmc.getX();

            double gma1 = gma.getY();
            double gmb1 = gmb.getY();
            double gmc1 = gmc.getY();

            double gma2 = gma.getZ();
            double gmb2 = gmb.getZ();
            double gmc2 = gmc.getZ();


            // get s-s-s element
            double expo = _decay_alpha * _decay_gamma * (_pos_alpha - _pos_gamma)*(_pos_alpha - _pos_gamma)
                    + _decay_gamma * _decay_gw * (_pos_gamma - _pos_gw) * (_pos_gamma - _pos_gw)
                    + _decay_alpha * _decay_gw * (_pos_alpha - _pos_gw) * (_pos_alpha - _pos_gw);

            const double pi = boost::math::constants::pi<double>();
            double prefak = pow(8.0 * _decay_alpha * _decay_gamma * _decay_gw / pi, 0.75) * pow(fak2, 1.5);

            double value = prefak * exp(-fak2 * expo);

            // check if it contributes
            if (value < gwaccuracy) return false;

            // cout << "             s-s-s: " << value << endl;
            //cout << "             I need a 3d array of sizes [" << _nalpha << ":" << _ngw << ":" << _ngamma << "]" << endl;
            // if it does, go on and create multiarray
            typedef boost::multi_array<double, 3> ma_type;
            typedef boost::multi_array_types::extent_range range;
            //typedef ma_type::index index;
            ma_type::extent_gen extents;
            ma_type S;
            S.resize(extents[ range(1, _nalpha + 1) ][ range(1, _ngw + 1) ][ range(1, _ngamma + 1)]);
            //cout << "                  and I got it!" << endl;
            //cout << "                  S-dim1 base:shape   [" << S.index_bases()[0] << ":" << S.shape()[0] << "]" << endl;
            //cout << "                  S-dim2 base:shape   [" << S.index_bases()[1] << ":" << S.shape()[1] << "]" << endl;
            //cout << "                  S-dim3 base:shape   [" << S.index_bases()[2] << ":" << S.shape()[2] << "]" << endl;



            // now fill
            S[1][1][1] = value;
             // cout << "                 set s-s-s value" << endl;
            if (_lmax_alpha >= 0 && _lmax_gw >= 1 && _lmax_gamma >= 0) {
                S[1][2][1] = gmc0 * S[1][1][1];
                S[1][3][1] = gmc1 * S[1][1][1];
                S[1][4][1] = gmc2 * S[1][1][1];
            }

            if (_lmax_alpha >= 0 && _lmax_gw >= 2 && _lmax_gamma >= 0) {
                S[1][8][1] = gmc0 * S[1][2][1] + fak * S[1][1][1];
                S[1][5][1] = gmc1 * S[1][2][1];
                S[1][6][1] = gmc2 * S[1][2][1];
                S[1][9][1] = gmc1 * S[1][3][1] + fak * S[1][1][1];
                S[1][7][1] = gmc2 * S[1][3][1];
                S[1][10][1] = gmc2 * S[1][4][1] + fak * S[1][1][1];
            }

            if (_lmax_alpha >= 1 && _lmax_gw >= 0 && _lmax_gamma >= 0) {
                S[2][1][1] = gma0 * S[1][1][1];
                S[3][1][1] = gma1 * S[1][1][1];
                S[4][1][1] = gma2 * S[1][1][1];
            }

            if (_lmax_alpha >= 1 && _lmax_gw >= 1 && _lmax_gamma >= 0) {
                S[2][2][1] = gma0 * S[1][2][1] + fak * S[1][1][1];
                S[3][2][1] = gma1 * S[1][2][1];
                S[4][2][1] = gma2 * S[1][2][1];
                S[2][3][1] = gma0 * S[1][3][1];
                S[3][3][1] = gma1 * S[1][3][1] + fak * S[1][1][1];
                S[4][3][1] = gma2 * S[1][3][1];
                S[2][4][1] = gma0 * S[1][4][1];
                S[3][4][1] = gma1 * S[1][4][1];
                S[4][4][1] = gma2 * S[1][4][1] + fak * S[1][1][1];
            }

            if (_lmax_alpha >= 1 && _lmax_gw >= 2 && _lmax_gamma >= 0) {
                S[2][5][1] = gma0 * S[1][5][1] + fak * S[1][3][1];
                S[3][5][1] = gma1 * S[1][5][1] + fak * S[1][2][1];
                S[4][5][1] = gma2 * S[1][5][1];
                S[2][6][1] = gma0 * S[1][6][1] + fak * S[1][4][1];
                S[3][6][1] = gma1 * S[1][6][1];
                S[4][6][1] = gma2 * S[1][6][1] + fak * S[1][2][1];
                S[2][7][1] = gma0 * S[1][7][1];
                S[3][7][1] = gma1 * S[1][7][1] + fak * S[1][4][1];
                S[4][7][1] = gma2 * S[1][7][1] + fak * S[1][3][1];
                S[2][8][1] = gma0 * S[1][8][1] + fak2 * S[1][2][1];
                S[3][8][1] = gma1 * S[1][8][1];
                S[4][8][1] = gma2 * S[1][8][1];
                S[2][9][1] = gma0 * S[1][9][1];
                S[3][9][1] = gma1 * S[1][9][1] + fak2 * S[1][3][1];
                S[4][9][1] = gma2 * S[1][9][1];
                S[2][10][1] = gma0 * S[1][10][1];
                S[3][10][1] = gma1 * S[1][10][1];
                S[4][10][1] = gma2 * S[1][10][1] + fak2 * S[1][4][1];
            }

            if (_lmax_alpha >= 0 && _lmax_gw >= 0 && _lmax_gamma >= 1) {
                S[1][1][2] = gmb0 * S[1][1][1];
                S[1][1][3] = gmb1 * S[1][1][1];
                S[1][1][4] = gmb2 * S[1][1][1];
            }

            if (_lmax_alpha >= 0 && _lmax_gw >= 1 && _lmax_gamma >= 1) {
                S[1][2][2] = gmc0 * S[1][1][2] + fak * S[1][1][1];
                S[1][3][2] = gmc1 * S[1][1][2];
                S[1][4][2] = gmc2 * S[1][1][2];
                S[1][2][3] = gmc0 * S[1][1][3];
                S[1][3][3] = gmc1 * S[1][1][3] + fak * S[1][1][1];
                S[1][4][3] = gmc2 * S[1][1][3];
                S[1][2][4] = gmc0 * S[1][1][4];
                S[1][3][4] = gmc1 * S[1][1][4];
                S[1][4][4] = gmc2 * S[1][1][4] + fak * S[1][1][1];
            }

            if (_lmax_alpha >= 0 && _lmax_gw >= 2 && _lmax_gamma >= 1) {
                S[1][8][2] = gmc0 * S[1][2][2] + fak * (S[1][1][2] + S[1][2][1]);
                S[1][5][2] = gmc1 * S[1][2][2];
                S[1][6][2] = gmc2 * S[1][2][2];
                S[1][8][3] = gmc0 * S[1][2][3] + fak * S[1][1][3];
                S[1][5][3] = gmc1 * S[1][2][3] + fak * S[1][2][1];
                S[1][6][3] = gmc2 * S[1][2][3];
                S[1][8][4] = gmc0 * S[1][2][4] + fak * S[1][1][4];
                S[1][5][4] = gmc1 * S[1][2][4];
                S[1][6][4] = gmc2 * S[1][2][4] + fak * S[1][2][1];
                S[1][9][2] = gmc1 * S[1][3][2] + fak * S[1][1][2];
                S[1][7][2] = gmc2 * S[1][3][2];
                S[1][9][3] = gmc1 * S[1][3][3] + fak * (S[1][1][3] + S[1][3][1]);
                S[1][7][3] = gmc2 * S[1][3][3];
                S[1][9][4] = gmc1 * S[1][3][4] + fak * S[1][1][4];
                S[1][7][4] = gmc2 * S[1][3][4] + fak * S[1][3][1];
                S[1][10][2] = gmc2 * S[1][4][2] + fak * S[1][1][2];
                S[1][10][3] = gmc2 * S[1][4][3] + fak * S[1][1][3];
                S[1][10][4] = gmc2 * S[1][4][4] + fak * (S[1][1][4] + S[1][4][1]);
            }

            if (_lmax_alpha >= 1 && _lmax_gw >= 0 && _lmax_gamma >= 1) {
                S[2][1][2] = gma0 * S[1][1][2] + fak * S[1][1][1];
                S[3][1][2] = gma1 * S[1][1][2];
                S[4][1][2] = gma2 * S[1][1][2];
                S[2][1][3] = gma0 * S[1][1][3];
                S[3][1][3] = gma1 * S[1][1][3] + fak * S[1][1][1];
                S[4][1][3] = gma2 * S[1][1][3];
                S[2][1][4] = gma0 * S[1][1][4];
                S[3][1][4] = gma1 * S[1][1][4];
                S[4][1][4] = gma2 * S[1][1][4] + fak * S[1][1][1];
            }

            if (_lmax_alpha >= 1 && _lmax_gw >= 1 && _lmax_gamma >= 1) {
                S[2][2][2] = gma0 * S[1][2][2] + fak * (S[1][1][2] + S[1][2][1]);
                S[3][2][2] = gma1 * S[1][2][2];
                S[4][2][2] = gma2 * S[1][2][2];
                S[2][2][3] = gma0 * S[1][2][3] + fak * S[1][1][3];
                S[3][2][3] = gma1 * S[1][2][3] + fak * S[1][2][1];
                S[4][2][3] = gma2 * S[1][2][3];
                S[2][2][4] = gma0 * S[1][2][4] + fak * S[1][1][4];
                S[3][2][4] = gma1 * S[1][2][4];
                S[4][2][4] = gma2 * S[1][2][4] + fak * S[1][2][1];
                S[2][3][2] = gma0 * S[1][3][2] + fak * S[1][3][1];
                S[3][3][2] = gma1 * S[1][3][2] + fak * S[1][1][2];
                S[4][3][2] = gma2 * S[1][3][2];
                S[2][3][3] = gma0 * S[1][3][3];
                S[3][3][3] = gma1 * S[1][3][3] + fak * (S[1][1][3] + S[1][3][1]);
                S[4][3][3] = gma2 * S[1][3][3];
                S[2][3][4] = gma0 * S[1][3][4];
                S[3][3][4] = gma1 * S[1][3][4] + fak * S[1][1][4];
                S[4][3][4] = gma2 * S[1][3][4] + fak * S[1][3][1];
                S[2][4][2] = gma0 * S[1][4][2] + fak * S[1][4][1];
                S[3][4][2] = gma1 * S[1][4][2];
                S[4][4][2] = gma2 * S[1][4][2] + fak * S[1][1][2];
                S[2][4][3] = gma0 * S[1][4][3];
                S[3][4][3] = gma1 * S[1][4][3] + fak * S[1][4][1];
                S[4][4][3] = gma2 * S[1][4][3] + fak * S[1][1][3];
                S[2][4][4] = gma0 * S[1][4][4];
                S[3][4][4] = gma1 * S[1][4][4];
                S[4][4][4] = gma2 * S[1][4][4] + fak * (S[1][1][4] + S[1][4][1]);
            }

            if (_lmax_alpha >= 1 && _lmax_gw >= 2 && _lmax_gamma >= 1) {
                S[2][5][2] = gma0 * S[1][5][2] + fak * (S[1][3][2] + S[1][5][1]);
                S[3][5][2] = gma1 * S[1][5][2] + fak * S[1][2][2];
                S[4][5][2] = gma2 * S[1][5][2];
                S[2][5][3] = gma0 * S[1][5][3] + fak * S[1][3][3];
                S[3][5][3] = gma1 * S[1][5][3] + fak * (S[1][2][3] + S[1][5][1]);
                S[4][5][3] = gma2 * S[1][5][3];
                S[2][5][4] = gma0 * S[1][5][4] + fak * S[1][3][4];
                S[3][5][4] = gma1 * S[1][5][4] + fak * S[1][2][4];
                S[4][5][4] = gma2 * S[1][5][4] + fak * S[1][5][1];
                S[2][6][2] = gma0 * S[1][6][2] + fak * (S[1][4][2] + S[1][6][1]);
                S[3][6][2] = gma1 * S[1][6][2];
                S[4][6][2] = gma2 * S[1][6][2] + fak * S[1][2][2];
                S[2][6][3] = gma0 * S[1][6][3] + fak * S[1][4][3];
                S[3][6][3] = gma1 * S[1][6][3] + fak * S[1][6][1];
                S[4][6][3] = gma2 * S[1][6][3] + fak * S[1][2][3];
                S[2][6][4] = gma0 * S[1][6][4] + fak * S[1][4][4];
                S[3][6][4] = gma1 * S[1][6][4];
                S[4][6][4] = gma2 * S[1][6][4] + fak * (S[1][2][4] + S[1][6][1]);
                S[2][7][2] = gma0 * S[1][7][2] + fak * S[1][7][1];
                S[3][7][2] = gma1 * S[1][7][2] + fak * S[1][4][2];
                S[4][7][2] = gma2 * S[1][7][2] + fak * S[1][3][2];
                S[2][7][3] = gma0 * S[1][7][3];
                S[3][7][3] = gma1 * S[1][7][3] + fak * (S[1][4][3] + S[1][7][1]);
                S[4][7][3] = gma2 * S[1][7][3] + fak * S[1][3][3];
                S[2][7][4] = gma0 * S[1][7][4];
                S[3][7][4] = gma1 * S[1][7][4] + fak * S[1][4][4];
                S[4][7][4] = gma2 * S[1][7][4] + fak * (S[1][3][4] + S[1][7][1]);
                S[2][8][2] = gma0 * S[1][8][2] + fak * (S[1][2][2] + S[1][2][2] + S[1][8][1]);
                S[3][8][2] = gma1 * S[1][8][2];
                S[4][8][2] = gma2 * S[1][8][2];
                S[2][8][3] = gma0 * S[1][8][3] + fak2 * S[1][2][3];
                S[3][8][3] = gma1 * S[1][8][3] + fak * S[1][8][1];
                S[4][8][3] = gma2 * S[1][8][3];
                S[2][8][4] = gma0 * S[1][8][4] + fak * (2 * S[1][2][4]);
                S[3][8][4] = gma1 * S[1][8][4];
                S[4][8][4] = gma2 * S[1][8][4] + fak * S[1][8][1];
                S[2][9][2] = gma0 * S[1][9][2] + fak * S[1][9][1];
                S[3][9][2] = gma1 * S[1][9][2] + fak2 * S[1][3][2];
                S[4][9][2] = gma2 * S[1][9][2];
                S[2][9][3] = gma0 * S[1][9][3];
                S[3][9][3] = gma1 * S[1][9][3] + fak * (S[1][3][3] + S[1][3][3] + S[1][9][1]);
                S[4][9][3] = gma2 * S[1][9][3];
                S[2][9][4] = gma0 * S[1][9][4];
                S[3][9][4] = gma1 * S[1][9][4] + fak2 * S[1][3][4];
                S[4][9][4] = gma2 * S[1][9][4] + fak * S[1][9][1];
                S[2][10][2] = gma0 * S[1][10][2] + fak * S[1][10][1];
                S[3][10][2] = gma1 * S[1][10][2];
                S[4][10][2] = gma2 * S[1][10][2] + fak2 * S[1][4][2];
                S[2][10][3] = gma0 * S[1][10][3];
                S[3][10][3] = gma1 * S[1][10][3] + fak * S[1][10][1];
                S[4][10][3] = gma2 * S[1][10][3] + fak2 * S[1][4][3];
                S[2][10][4] = gma0 * S[1][10][4];
                S[3][10][4] = gma1 * S[1][10][4];
                S[4][10][4] = gma2 * S[1][10][4] + fak * (S[1][4][4] + S[1][4][4] + S[1][10][1]);
            }

            if (_lmax_alpha >= 0 && _lmax_gw >= 0 && _lmax_gamma >= 2) {
                S[1][1][8] = gmb0 * S[1][1][2] + fak * S[1][1][1];
                S[1][1][5] = gmb1 * S[1][1][2];
                S[1][1][6] = gmb2 * S[1][1][2];
                S[1][1][9] = gmb1 * S[1][1][3] + fak * S[1][1][1];
                S[1][1][7] = gmb2 * S[1][1][3];
                S[1][1][10] = gmb2 * S[1][1][4] + fak * S[1][1][1];
            }

            if (_lmax_alpha >= 0 && _lmax_gw >= 1 && _lmax_gamma >= 2) {
                S[1][2][5] = gmc0 * S[1][1][5] + fak * S[1][1][3];
                S[1][3][5] = gmc1 * S[1][1][5] + fak * S[1][1][2];
                S[1][4][5] = gmc2 * S[1][1][5];
                S[1][2][6] = gmc0 * S[1][1][6] + fak * S[1][1][4];
                S[1][3][6] = gmc1 * S[1][1][6];
                S[1][4][6] = gmc2 * S[1][1][6] + fak * S[1][1][2];
                S[1][2][7] = gmc0 * S[1][1][7];
                S[1][3][7] = gmc1 * S[1][1][7] + fak * S[1][1][4];
                S[1][4][7] = gmc2 * S[1][1][7] + fak * S[1][1][3];
                S[1][2][8] = gmc0 * S[1][1][8] + fak2 * S[1][1][2];
                S[1][3][8] = gmc1 * S[1][1][8];
                S[1][4][8] = gmc2 * S[1][1][8];
                S[1][2][9] = gmc0 * S[1][1][9];
                S[1][3][9] = gmc1 * S[1][1][9] + fak2 * S[1][1][3];
                S[1][4][9] = gmc2 * S[1][1][9];
                S[1][2][10] = gmc0 * S[1][1][10];
                S[1][3][10] = gmc1 * S[1][1][10];
                S[1][4][10] = gmc2 * S[1][1][10] + fak2 * S[1][1][4];
            }

            if (_lmax_alpha >= 0 && _lmax_gw >= 2 && _lmax_gamma >= 2) {
                S[1][8][5] = gmc0 * S[1][2][5] + fak * (S[1][1][5] + S[1][2][3]);
                S[1][5][5] = gmc1 * S[1][2][5] + fak * S[1][2][2];
                S[1][6][5] = gmc2 * S[1][2][5];
                S[1][8][6] = gmc0 * S[1][2][6] + fak * (S[1][1][6] + S[1][2][4]);
                S[1][5][6] = gmc1 * S[1][2][6];
                S[1][6][6] = gmc2 * S[1][2][6] + fak * S[1][2][2];
                S[1][8][7] = gmc0 * S[1][2][7] + fak * S[1][1][7];
                S[1][5][7] = gmc1 * S[1][2][7] + fak * S[1][2][4];
                S[1][6][7] = gmc2 * S[1][2][7] + fak * S[1][2][3];
                S[1][8][8] = gmc0 * S[1][2][8] + fak * (S[1][1][8] + S[1][2][2] + S[1][2][2]);
                S[1][5][8] = gmc1 * S[1][2][8];
                S[1][6][8] = gmc2 * S[1][2][8];
                S[1][8][9] = gmc0 * S[1][2][9] + fak * S[1][1][9];
                S[1][5][9] = gmc1 * S[1][2][9] + fak2 * S[1][2][3];
                S[1][6][9] = gmc2 * S[1][2][9];
                S[1][8][10] = gmc0 * S[1][2][10] + fak * S[1][1][10];
                S[1][5][10] = gmc1 * S[1][2][10];
                S[1][6][10] = gmc2 * S[1][2][10] + fak2 * S[1][2][4];
                S[1][9][5] = gmc1 * S[1][3][5] + fak * (S[1][1][5] + S[1][3][2]);
                S[1][7][5] = gmc2 * S[1][3][5];
                S[1][9][6] = gmc1 * S[1][3][6] + fak * S[1][1][6];
                S[1][7][6] = gmc2 * S[1][3][6] + fak * S[1][3][2];
                S[1][9][7] = gmc1 * S[1][3][7] + fak * (S[1][1][7] + S[1][3][4]);
                S[1][7][7] = gmc2 * S[1][3][7] + fak * S[1][3][3];
                S[1][9][8] = gmc1 * S[1][3][8] + fak * S[1][1][8];
                S[1][7][8] = gmc2 * S[1][3][8];
                S[1][9][9] = gmc1 * S[1][3][9] + fak * (S[1][1][9] + S[1][3][3] + S[1][3][3]);
                S[1][7][9] = gmc2 * S[1][3][9];
                S[1][9][10] = gmc1 * S[1][3][10] + fak * S[1][1][10];
                S[1][7][10] = gmc2 * S[1][3][10] + fak2 * S[1][3][4];
                S[1][10][5] = gmc2 * S[1][4][5] + fak * S[1][1][5];
                S[1][10][6] = gmc2 * S[1][4][6] + fak * (S[1][1][6] + S[1][4][2]);
                S[1][10][7] = gmc2 * S[1][4][7] + fak * (S[1][1][7] + S[1][4][3]);
                S[1][10][8] = gmc2 * S[1][4][8] + fak * S[1][1][8];
                S[1][10][9] = gmc2 * S[1][4][9] + fak * S[1][1][9];
                S[1][10][10] = gmc2 * S[1][4][10] + fak * (S[1][1][10] + S[1][4][4] + S[1][4][4]);
            }

            if (_lmax_alpha >= 2 && _lmax_gw >= 0 && _lmax_gamma >= 0) {
                S[8][1][1] = gma0 * S[2][1][1] + fak * S[1][1][1];
                S[5][1][1] = gma1 * S[2][1][1];
                S[6][1][1] = gma2 * S[2][1][1];
                S[9][1][1] = gma1 * S[3][1][1] + fak * S[1][1][1];
                S[7][1][1] = gma2 * S[3][1][1];
                S[10][1][1] = gma2 * S[4][1][1] + fak * S[1][1][1];
            }

            if (_lmax_alpha >= 2 && _lmax_gw >= 1 && _lmax_gamma >= 0) {
                S[8][2][1] = gma0 * S[2][2][1] + fak * (S[1][2][1] + S[2][1][1]);
                S[5][2][1] = gma1 * S[2][2][1];
                S[6][2][1] = gma2 * S[2][2][1];
                S[8][3][1] = gma0 * S[2][3][1] + fak * S[1][3][1];
                S[5][3][1] = gma1 * S[2][3][1] + fak * S[2][1][1];
                S[6][3][1] = gma2 * S[2][3][1];
                S[8][4][1] = gma0 * S[2][4][1] + fak * S[1][4][1];
                S[5][4][1] = gma1 * S[2][4][1];
                S[6][4][1] = gma2 * S[2][4][1] + fak * S[2][1][1];
                S[9][2][1] = gma1 * S[3][2][1] + fak * S[1][2][1];
                S[7][2][1] = gma2 * S[3][2][1];
                S[9][3][1] = gma1 * S[3][3][1] + fak * (S[1][3][1] + S[3][1][1]);
                S[7][3][1] = gma2 * S[3][3][1];
                S[9][4][1] = gma1 * S[3][4][1] + fak * S[1][4][1];
                S[7][4][1] = gma2 * S[3][4][1] + fak * S[3][1][1];
                S[10][2][1] = gma2 * S[4][2][1] + fak * S[1][2][1];
                S[10][3][1] = gma2 * S[4][3][1] + fak * S[1][3][1];
                S[10][4][1] = gma2 * S[4][4][1] + fak * (S[1][4][1] + S[4][1][1]);
            }

            if (_lmax_alpha >= 2 && _lmax_gw >= 2 && _lmax_gamma >= 0) {
                S[8][5][1] = gma0 * S[2][5][1] + fak * (S[1][5][1] + S[2][3][1]);
                S[5][5][1] = gma1 * S[2][5][1] + fak * S[2][2][1];
                S[6][5][1] = gma2 * S[2][5][1];
                S[8][6][1] = gma0 * S[2][6][1] + fak * (S[1][6][1] + S[2][4][1]);
                S[5][6][1] = gma1 * S[2][6][1];
                S[6][6][1] = gma2 * S[2][6][1] + fak * S[2][2][1];
                S[8][7][1] = gma0 * S[2][7][1] + fak * S[1][7][1];
                S[5][7][1] = gma1 * S[2][7][1] + fak * S[2][4][1];
                S[6][7][1] = gma2 * S[2][7][1] + fak * S[2][3][1];
                S[8][8][1] = gma0 * S[2][8][1] + fak * (S[1][8][1] + S[2][2][1] + S[2][2][1]);
                S[5][8][1] = gma1 * S[2][8][1];
                S[6][8][1] = gma2 * S[2][8][1];
                S[8][9][1] = gma0 * S[2][9][1] + fak * S[1][9][1];
                S[5][9][1] = gma1 * S[2][9][1] + fak2 * S[2][3][1];
                S[6][9][1] = gma2 * S[2][9][1];
                S[8][10][1] = gma0 * S[2][10][1] + fak * S[1][10][1];
                S[5][10][1] = gma1 * S[2][10][1];
                S[6][10][1] = gma2 * S[2][10][1] + fak2 * S[2][4][1];
                S[9][5][1] = gma1 * S[3][5][1] + fak * (S[1][5][1] + S[3][2][1]);
                S[7][5][1] = gma2 * S[3][5][1];
                S[9][6][1] = gma1 * S[3][6][1] + fak * S[1][6][1];
                S[7][6][1] = gma2 * S[3][6][1] + fak * S[3][2][1];
                S[9][7][1] = gma1 * S[3][7][1] + fak * (S[1][7][1] + S[3][4][1]);
                S[7][7][1] = gma2 * S[3][7][1] + fak * S[3][3][1];
                S[9][8][1] = gma1 * S[3][8][1] + fak * S[1][8][1];
                S[7][8][1] = gma2 * S[3][8][1];
                S[9][9][1] = gma1 * S[3][9][1] + fak * (S[1][9][1] + S[3][3][1] + S[3][3][1]);
                S[7][9][1] = gma2 * S[3][9][1];
                S[9][10][1] = gma1 * S[3][10][1] + fak * S[1][10][1];
                S[7][10][1] = gma2 * S[3][10][1] + fak2 * S[3][4][1];
                S[10][5][1] = gma2 * S[4][5][1] + fak * S[1][5][1];
                S[10][6][1] = gma2 * S[4][6][1] + fak * (S[1][6][1] + S[4][2][1]);
                S[10][7][1] = gma2 * S[4][7][1] + fak * (S[1][7][1] + S[4][3][1]);
                S[10][8][1] = gma2 * S[4][8][1] + fak * S[1][8][1];
                S[10][9][1] = gma2 * S[4][9][1] + fak * S[1][9][1];
                S[10][10][1] = gma2 * S[4][10][1] + fak * (S[1][10][1] + S[4][4][1] + S[4][4][1]);
            }

            if (_lmax_alpha >= 1 && _lmax_gw >= 0 && _lmax_gamma >= 2) {
                S[2][1][5] = gma0 * S[1][1][5] + fak * S[1][1][3];
                S[3][1][5] = gma1 * S[1][1][5] + fak * S[1][1][2];
                S[4][1][5] = gma2 * S[1][1][5];
                S[2][1][6] = gma0 * S[1][1][6] + fak * S[1][1][4];
                S[3][1][6] = gma1 * S[1][1][6];
                S[4][1][6] = gma2 * S[1][1][6] + fak * S[1][1][2];
                S[2][1][7] = gma0 * S[1][1][7];
                S[3][1][7] = gma1 * S[1][1][7] + fak * S[1][1][4];
                S[4][1][7] = gma2 * S[1][1][7] + fak * S[1][1][3];
                S[2][1][8] = gma0 * S[1][1][8] + fak2 * S[1][1][2];
                S[3][1][8] = gma1 * S[1][1][8];
                S[4][1][8] = gma2 * S[1][1][8];
                S[2][1][9] = gma0 * S[1][1][9];
                S[3][1][9] = gma1 * S[1][1][9] + fak2 * S[1][1][3];
                S[4][1][9] = gma2 * S[1][1][9];
                S[2][1][10] = gma0 * S[1][1][10];
                S[3][1][10] = gma1 * S[1][1][10];
                S[4][1][10] = gma2 * S[1][1][10] + fak2 * S[1][1][4];
            }

            if (_lmax_alpha >= 1 && _lmax_gw >= 1 && _lmax_gamma >= 2) {
                S[2][2][5] = gma0 * S[1][2][5] + fak * (S[1][1][5] + S[1][2][3]);
                S[3][2][5] = gma1 * S[1][2][5] + fak * S[1][2][2];
                S[4][2][5] = gma2 * S[1][2][5];
                S[2][2][6] = gma0 * S[1][2][6] + fak * (S[1][1][6] + S[1][2][4]);
                S[3][2][6] = gma1 * S[1][2][6];
                S[4][2][6] = gma2 * S[1][2][6] + fak * S[1][2][2];
                S[2][2][7] = gma0 * S[1][2][7] + fak * S[1][1][7];
                S[3][2][7] = gma1 * S[1][2][7] + fak * S[1][2][4];
                S[4][2][7] = gma2 * S[1][2][7] + fak * S[1][2][3];
                S[2][2][8] = gma0 * S[1][2][8] + fak * (S[1][1][8] + S[1][2][2] + S[1][2][2]);
                S[3][2][8] = gma1 * S[1][2][8];
                S[4][2][8] = gma2 * S[1][2][8];
                S[2][2][9] = gma0 * S[1][2][9] + fak * S[1][1][9];
                S[3][2][9] = gma1 * S[1][2][9] + fak2 * S[1][2][3];
                S[4][2][9] = gma2 * S[1][2][9];
                S[2][2][10] = gma0 * S[1][2][10] + fak * S[1][1][10];
                S[3][2][10] = gma1 * S[1][2][10];
                S[4][2][10] = gma2 * S[1][2][10] + fak2 * S[1][2][4];
                S[2][3][5] = gma0 * S[1][3][5] + fak * S[1][3][3];
                S[3][3][5] = gma1 * S[1][3][5] + fak * (S[1][1][5] + S[1][3][2]);
                S[4][3][5] = gma2 * S[1][3][5];
                S[2][3][6] = gma0 * S[1][3][6] + fak * S[1][3][4];
                S[3][3][6] = gma1 * S[1][3][6] + fak * S[1][1][6];
                S[4][3][6] = gma2 * S[1][3][6] + fak * S[1][3][2];
                S[2][3][7] = gma0 * S[1][3][7];
                S[3][3][7] = gma1 * S[1][3][7] + fak * (S[1][1][7] + S[1][3][4]);
                S[4][3][7] = gma2 * S[1][3][7] + fak * S[1][3][3];
                S[2][3][8] = gma0 * S[1][3][8] + fak2 * S[1][3][2];
                S[3][3][8] = gma1 * S[1][3][8] + fak * S[1][1][8];
                S[4][3][8] = gma2 * S[1][3][8];
                S[2][3][9] = gma0 * S[1][3][9];
                S[3][3][9] = gma1 * S[1][3][9] + fak * (S[1][1][9] + S[1][3][3] + S[1][3][3]);
                S[4][3][9] = gma2 * S[1][3][9];
                S[2][3][10] = gma0 * S[1][3][10];
                S[3][3][10] = gma1 * S[1][3][10] + fak * S[1][1][10];
                S[4][3][10] = gma2 * S[1][3][10] + fak2 * S[1][3][4];
                S[2][4][5] = gma0 * S[1][4][5] + fak * S[1][4][3];
                S[3][4][5] = gma1 * S[1][4][5] + fak * S[1][4][2];
                S[4][4][5] = gma2 * S[1][4][5] + fak * S[1][1][5];
                S[2][4][6] = gma0 * S[1][4][6] + fak * S[1][4][4];
                S[3][4][6] = gma1 * S[1][4][6];
                S[4][4][6] = gma2 * S[1][4][6] + fak * (S[1][1][6] + S[1][4][2]);
                S[2][4][7] = gma0 * S[1][4][7];
                S[3][4][7] = gma1 * S[1][4][7] + fak * S[1][4][4];
                S[4][4][7] = gma2 * S[1][4][7] + fak * (S[1][1][7] + S[1][4][3]);
                S[2][4][8] = gma0 * S[1][4][8] + fak2 * S[1][4][2];
                S[3][4][8] = gma1 * S[1][4][8];
                S[4][4][8] = gma2 * S[1][4][8] + fak * S[1][1][8];
                S[2][4][9] = gma0 * S[1][4][9];
                S[3][4][9] = gma1 * S[1][4][9] + fak2 * S[1][4][3];
                S[4][4][9] = gma2 * S[1][4][9] + fak * S[1][1][9];
                S[2][4][10] = gma0 * S[1][4][10];
                S[3][4][10] = gma1 * S[1][4][10];
                S[4][4][10] = gma2 * S[1][4][10] + fak * (S[1][1][10] + S[1][4][4] + S[1][4][4]);
            }

            if (_lmax_alpha >= 1 && _lmax_gw >= 2 && _lmax_gamma >= 2) {
                S[2][5][5] = gma0 * S[1][5][5] + fak * (S[1][3][5] + S[1][5][3]);
                S[3][5][5] = gma1 * S[1][5][5] + fak * (S[1][2][5] + S[1][5][2]);
                S[4][5][5] = gma2 * S[1][5][5];
                S[2][5][6] = gma0 * S[1][5][6] + fak * (S[1][3][6] + S[1][5][4]);
                S[3][5][6] = gma1 * S[1][5][6] + fak * S[1][2][6];
                S[4][5][6] = gma2 * S[1][5][6] + fak * S[1][5][2];
                S[2][5][7] = gma0 * S[1][5][7] + fak * S[1][3][7];
                S[3][5][7] = gma1 * S[1][5][7] + fak * (S[1][2][7] + S[1][5][4]);
                S[4][5][7] = gma2 * S[1][5][7] + fak * S[1][5][3];
                S[2][5][8] = gma0 * S[1][5][8] + fak * (S[1][3][8] + S[1][5][2] + S[1][5][2]);
                S[3][5][8] = gma1 * S[1][5][8] + fak * S[1][2][8];
                S[4][5][8] = gma2 * S[1][5][8];
                S[2][5][9] = gma0 * S[1][5][9] + fak * S[1][3][9];
                S[3][5][9] = gma1 * S[1][5][9] + fak * (S[1][2][9] + S[1][5][3] + S[1][5][3]);
                S[4][5][9] = gma2 * S[1][5][9];
                S[2][5][10] = gma0 * S[1][5][10] + fak * S[1][3][10];
                S[3][5][10] = gma1 * S[1][5][10] + fak * S[1][2][10];
                S[4][5][10] = gma2 * S[1][5][10] + fak2 * S[1][5][4];
                S[2][6][5] = gma0 * S[1][6][5] + fak * (S[1][4][5] + S[1][6][3]);
                S[3][6][5] = gma1 * S[1][6][5] + fak * S[1][6][2];
                S[4][6][5] = gma2 * S[1][6][5] + fak * S[1][2][5];
                S[2][6][6] = gma0 * S[1][6][6] + fak * (S[1][4][6] + S[1][6][4]);
                S[3][6][6] = gma1 * S[1][6][6];
                S[4][6][6] = gma2 * S[1][6][6] + fak * (S[1][2][6] + S[1][6][2]);
                S[2][6][7] = gma0 * S[1][6][7] + fak * S[1][4][7];
                S[3][6][7] = gma1 * S[1][6][7] + fak * S[1][6][4];
                S[4][6][7] = gma2 * S[1][6][7] + fak * (S[1][2][7] + S[1][6][3]);
                S[2][6][8] = gma0 * S[1][6][8] + fak * (S[1][4][8] + S[1][6][2] + S[1][6][2]);
                S[3][6][8] = gma1 * S[1][6][8];
                S[4][6][8] = gma2 * S[1][6][8] + fak * S[1][2][8];
                S[2][6][9] = gma0 * S[1][6][9] + fak * S[1][4][9];
                S[3][6][9] = gma1 * S[1][6][9] + fak2 * S[1][6][3];
                S[4][6][9] = gma2 * S[1][6][9] + fak * S[1][2][9];
                S[2][6][10] = gma0 * S[1][6][10] + fak * S[1][4][10];
                S[3][6][10] = gma1 * S[1][6][10];
                S[4][6][10] = gma2 * S[1][6][10] + fak * (S[1][2][10] + S[1][6][4] + S[1][6][4]);
                S[2][7][5] = gma0 * S[1][7][5] + fak * S[1][7][3];
                S[3][7][5] = gma1 * S[1][7][5] + fak * (S[1][4][5] + S[1][7][2]);
                S[4][7][5] = gma2 * S[1][7][5] + fak * S[1][3][5];
                S[2][7][6] = gma0 * S[1][7][6] + fak * S[1][7][4];
                S[3][7][6] = gma1 * S[1][7][6] + fak * S[1][4][6];
                S[4][7][6] = gma2 * S[1][7][6] + fak * (S[1][3][6] + S[1][7][2]);
                S[2][7][7] = gma0 * S[1][7][7];
                S[3][7][7] = gma1 * S[1][7][7] + fak * (S[1][4][7] + S[1][7][4]);
                S[4][7][7] = gma2 * S[1][7][7] + fak * (S[1][3][7] + S[1][7][3]);
                S[2][7][8] = gma0 * S[1][7][8] + fak * (2 * S[1][7][2]);
                S[3][7][8] = gma1 * S[1][7][8] + fak * S[1][4][8];
                S[4][7][8] = gma2 * S[1][7][8] + fak * S[1][3][8];
                S[2][7][9] = gma0 * S[1][7][9];
                S[3][7][9] = gma1 * S[1][7][9] + fak * (S[1][4][9] + S[1][7][3] + S[1][7][3]);
                S[4][7][9] = gma2 * S[1][7][9] + fak * S[1][3][9];
                S[2][7][10] = gma0 * S[1][7][10];
                S[3][7][10] = gma1 * S[1][7][10] + fak * S[1][4][10];
                S[4][7][10] = gma2 * S[1][7][10] + fak * (S[1][3][10] + S[1][7][4] + S[1][7][4]);
                S[2][8][5] = gma0 * S[1][8][5] + fak * (S[1][2][5] + S[1][2][5] + S[1][8][3]);
                S[3][8][5] = gma1 * S[1][8][5] + fak * S[1][8][2];
                S[4][8][5] = gma2 * S[1][8][5];
                S[2][8][6] = gma0 * S[1][8][6] + fak * (S[1][2][6] + S[1][2][6] + S[1][8][4]);
                S[3][8][6] = gma1 * S[1][8][6];
                S[4][8][6] = gma2 * S[1][8][6] + fak * S[1][8][2];
                S[2][8][7] = gma0 * S[1][8][7] + fak2 * S[1][2][7];
                S[3][8][7] = gma1 * S[1][8][7] + fak * S[1][8][4];
                S[4][8][7] = gma2 * S[1][8][7] + fak * S[1][8][3];
                S[2][8][8] = gma0 * S[1][8][8] + fak2 * (S[1][2][8] + S[1][8][2]);
                S[3][8][8] = gma1 * S[1][8][8];
                S[4][8][8] = gma2 * S[1][8][8];
                S[2][8][9] = gma0 * S[1][8][9] + fak2 * S[1][2][9];
                S[3][8][9] = gma1 * S[1][8][9] + fak2 * S[1][8][3];
                S[4][8][9] = gma2 * S[1][8][9];
                S[2][8][10] = gma0 * S[1][8][10] + fak2 * S[1][2][10];
                S[3][8][10] = gma1 * S[1][8][10];
                S[4][8][10] = gma2 * S[1][8][10] + fak2 * S[1][8][4];
                S[2][9][5] = gma0 * S[1][9][5] + fak * S[1][9][3];
                S[3][9][5] = gma1 * S[1][9][5] + fak * (S[1][3][5] + S[1][3][5] + S[1][9][2]);
                S[4][9][5] = gma2 * S[1][9][5];
                S[2][9][6] = gma0 * S[1][9][6] + fak * S[1][9][4];
                S[3][9][6] = gma1 * S[1][9][6] + fak2 * S[1][3][6];
                S[4][9][6] = gma2 * S[1][9][6] + fak * S[1][9][2];
                S[2][9][7] = gma0 * S[1][9][7];
                S[3][9][7] = gma1 * S[1][9][7] + fak * (S[1][3][7] + S[1][3][7] + S[1][9][4]);
                S[4][9][7] = gma2 * S[1][9][7] + fak * S[1][9][3];
                S[2][9][8] = gma0 * S[1][9][8] + fak2 * S[1][9][2];
                S[3][9][8] = gma1 * S[1][9][8] + fak2 * S[1][3][8];
                S[4][9][8] = gma2 * S[1][9][8];
                S[2][9][9] = gma0 * S[1][9][9];
                S[3][9][9] = gma1 * S[1][9][9] + fak2 * (S[1][3][9] + S[1][9][3]);
                S[4][9][9] = gma2 * S[1][9][9];
                S[2][9][10] = gma0 * S[1][9][10];
                S[3][9][10] = gma1 * S[1][9][10] + fak2 * S[1][3][10];
                S[4][9][10] = gma2 * S[1][9][10] + fak2 * S[1][9][4];
                S[2][10][5] = gma0 * S[1][10][5] + fak * S[1][10][3];
                S[3][10][5] = gma1 * S[1][10][5] + fak * S[1][10][2];
                S[4][10][5] = gma2 * S[1][10][5] + fak2 * S[1][4][5];
                S[2][10][6] = gma0 * S[1][10][6] + fak * S[1][10][4];
                S[3][10][6] = gma1 * S[1][10][6];
                S[4][10][6] = gma2 * S[1][10][6] + fak * (S[1][4][6] + S[1][4][6] + S[1][10][2]);
                S[2][10][7] = gma0 * S[1][10][7];
                S[3][10][7] = gma1 * S[1][10][7] + fak * S[1][10][4];
                S[4][10][7] = gma2 * S[1][10][7] + fak * (S[1][4][7] + S[1][4][7] + S[1][10][3]);
                S[2][10][8] = gma0 * S[1][10][8] + fak2 * S[1][10][2];
                S[3][10][8] = gma1 * S[1][10][8];
                S[4][10][8] = gma2 * S[1][10][8] + fak2 * S[1][4][8];
                S[2][10][9] = gma0 * S[1][10][9];
                S[3][10][9] = gma1 * S[1][10][9] + fak2 * S[1][10][3];
                S[4][10][9] = gma2 * S[1][10][9] + fak2 * S[1][4][9];
                S[2][10][10] = gma0 * S[1][10][10];
                S[3][10][10] = gma1 * S[1][10][10];
                S[4][10][10] = gma2 * S[1][10][10] + fak2 * (S[1][4][10] + S[1][10][4]);
            }

            if (_lmax_alpha >= 2 && _lmax_gw >= 0 && _lmax_gamma >= 1) {
                S[8][1][2] = gma0 * S[2][1][2] + fak * (S[1][1][2] + S[2][1][1]);
                S[5][1][2] = gma1 * S[2][1][2];
                S[6][1][2] = gma2 * S[2][1][2];
                S[8][1][3] = gma0 * S[2][1][3] + fak * S[1][1][3];
                S[5][1][3] = gma1 * S[2][1][3] + fak * S[2][1][1];
                S[6][1][3] = gma2 * S[2][1][3];
                S[8][1][4] = gma0 * S[2][1][4] + fak * S[1][1][4];
                S[5][1][4] = gma1 * S[2][1][4];
                S[6][1][4] = gma2 * S[2][1][4] + fak * S[2][1][1];
                S[9][1][2] = gma1 * S[3][1][2] + fak * S[1][1][2];
                S[7][1][2] = gma2 * S[3][1][2];
                S[9][1][3] = gma1 * S[3][1][3] + fak * (S[1][1][3] + S[3][1][1]);
                S[7][1][3] = gma2 * S[3][1][3];
                S[9][1][4] = gma1 * S[3][1][4] + fak * S[1][1][4];
                S[7][1][4] = gma2 * S[3][1][4] + fak * S[3][1][1];
                S[10][1][2] = gma2 * S[4][1][2] + fak * S[1][1][2];
                S[10][1][3] = gma2 * S[4][1][3] + fak * S[1][1][3];
                S[10][1][4] = gma2 * S[4][1][4] + fak * (S[1][1][4] + S[4][1][1]);
            }

            if (_lmax_alpha >= 2 && _lmax_gw >= 1 && _lmax_gamma >= 1) {
                S[8][2][2] = gma0 * S[2][2][2] + fak * (S[1][2][2] + S[2][1][2] + S[2][2][1]);
                S[5][2][2] = gma1 * S[2][2][2];
                S[6][2][2] = gma2 * S[2][2][2];
                S[8][2][3] = gma0 * S[2][2][3] + fak * (S[1][2][3] + S[2][1][3]);
                S[5][2][3] = gma1 * S[2][2][3] + fak * S[2][2][1];
                S[6][2][3] = gma2 * S[2][2][3];
                S[8][2][4] = gma0 * S[2][2][4] + fak * (S[1][2][4] + S[2][1][4]);
                S[5][2][4] = gma1 * S[2][2][4];
                S[6][2][4] = gma2 * S[2][2][4] + fak * S[2][2][1];
                S[8][3][2] = gma0 * S[2][3][2] + fak * (S[1][3][2] + S[2][3][1]);
                S[5][3][2] = gma1 * S[2][3][2] + fak * S[2][1][2];
                S[6][3][2] = gma2 * S[2][3][2];
                S[8][3][3] = gma0 * S[2][3][3] + fak * S[1][3][3];
                S[5][3][3] = gma1 * S[2][3][3] + fak * (S[2][1][3] + S[2][3][1]);
                S[6][3][3] = gma2 * S[2][3][3];
                S[8][3][4] = gma0 * S[2][3][4] + fak * S[1][3][4];
                S[5][3][4] = gma1 * S[2][3][4] + fak * S[2][1][4];
                S[6][3][4] = gma2 * S[2][3][4] + fak * S[2][3][1];
                S[8][4][2] = gma0 * S[2][4][2] + fak * (S[1][4][2] + S[2][4][1]);
                S[5][4][2] = gma1 * S[2][4][2];
                S[6][4][2] = gma2 * S[2][4][2] + fak * S[2][1][2];
                S[8][4][3] = gma0 * S[2][4][3] + fak * S[1][4][3];
                S[5][4][3] = gma1 * S[2][4][3] + fak * S[2][4][1];
                S[6][4][3] = gma2 * S[2][4][3] + fak * S[2][1][3];
                S[8][4][4] = gma0 * S[2][4][4] + fak * S[1][4][4];
                S[5][4][4] = gma1 * S[2][4][4];
                S[6][4][4] = gma2 * S[2][4][4] + fak * (S[2][1][4] + S[2][4][1]);
                S[9][2][2] = gma1 * S[3][2][2] + fak * S[1][2][2];
                S[7][2][2] = gma2 * S[3][2][2];
                S[9][2][3] = gma1 * S[3][2][3] + fak * (S[1][2][3] + S[3][2][1]);
                S[7][2][3] = gma2 * S[3][2][3];
                S[9][2][4] = gma1 * S[3][2][4] + fak * S[1][2][4];
                S[7][2][4] = gma2 * S[3][2][4] + fak * S[3][2][1];
                S[9][3][2] = gma1 * S[3][3][2] + fak * (S[1][3][2] + S[3][1][2]);
                S[7][3][2] = gma2 * S[3][3][2];
                S[9][3][3] = gma1 * S[3][3][3] + fak * (S[1][3][3] + S[3][1][3] + S[3][3][1]);
                S[7][3][3] = gma2 * S[3][3][3];
                S[9][3][4] = gma1 * S[3][3][4] + fak * (S[1][3][4] + S[3][1][4]);
                S[7][3][4] = gma2 * S[3][3][4] + fak * S[3][3][1];
                S[9][4][2] = gma1 * S[3][4][2] + fak * S[1][4][2];
                S[7][4][2] = gma2 * S[3][4][2] + fak * S[3][1][2];
                S[9][4][3] = gma1 * S[3][4][3] + fak * (S[1][4][3] + S[3][4][1]);
                S[7][4][3] = gma2 * S[3][4][3] + fak * S[3][1][3];
                S[9][4][4] = gma1 * S[3][4][4] + fak * S[1][4][4];
                S[7][4][4] = gma2 * S[3][4][4] + fak * (S[3][1][4] + S[3][4][1]);
                S[10][2][2] = gma2 * S[4][2][2] + fak * S[1][2][2];
                S[10][2][3] = gma2 * S[4][2][3] + fak * S[1][2][3];
                S[10][2][4] = gma2 * S[4][2][4] + fak * (S[1][2][4] + S[4][2][1]);
                S[10][3][2] = gma2 * S[4][3][2] + fak * S[1][3][2];
                S[10][3][3] = gma2 * S[4][3][3] + fak * S[1][3][3];
                S[10][3][4] = gma2 * S[4][3][4] + fak * (S[1][3][4] + S[4][3][1]);
                S[10][4][2] = gma2 * S[4][4][2] + fak * (S[1][4][2] + S[4][1][2]);
                S[10][4][3] = gma2 * S[4][4][3] + fak * (S[1][4][3] + S[4][1][3]);
                S[10][4][4] = gma2 * S[4][4][4] + fak * (S[1][4][4] + S[4][1][4] + S[4][4][1]);
            }

            if (_lmax_alpha >= 2 && _lmax_gw >= 2 && _lmax_gamma >= 1) {
                S[8][5][2] = gma0 * S[2][5][2] + fak * (S[1][5][2] + S[2][3][2] + S[2][5][1]);
                S[5][5][2] = gma1 * S[2][5][2] + fak * S[2][2][2];
                S[6][5][2] = gma2 * S[2][5][2];
                S[8][5][3] = gma0 * S[2][5][3] + fak * (S[1][5][3] + S[2][3][3]);
                S[5][5][3] = gma1 * S[2][5][3] + fak * (S[2][2][3] + S[2][5][1]);
                S[6][5][3] = gma2 * S[2][5][3];
                S[8][5][4] = gma0 * S[2][5][4] + fak * (S[1][5][4] + S[2][3][4]);
                S[5][5][4] = gma1 * S[2][5][4] + fak * S[2][2][4];
                S[6][5][4] = gma2 * S[2][5][4] + fak * S[2][5][1];
                S[8][6][2] = gma0 * S[2][6][2] + fak * (S[1][6][2] + S[2][4][2] + S[2][6][1]);
                S[5][6][2] = gma1 * S[2][6][2];
                S[6][6][2] = gma2 * S[2][6][2] + fak * S[2][2][2];
                S[8][6][3] = gma0 * S[2][6][3] + fak * (S[1][6][3] + S[2][4][3]);
                S[5][6][3] = gma1 * S[2][6][3] + fak * S[2][6][1];
                S[6][6][3] = gma2 * S[2][6][3] + fak * S[2][2][3];
                S[8][6][4] = gma0 * S[2][6][4] + fak * (S[1][6][4] + S[2][4][4]);
                S[5][6][4] = gma1 * S[2][6][4];
                S[6][6][4] = gma2 * S[2][6][4] + fak * (S[2][2][4] + S[2][6][1]);
                S[8][7][2] = gma0 * S[2][7][2] + fak * (S[1][7][2] + S[2][7][1]);
                S[5][7][2] = gma1 * S[2][7][2] + fak * S[2][4][2];
                S[6][7][2] = gma2 * S[2][7][2] + fak * S[2][3][2];
                S[8][7][3] = gma0 * S[2][7][3] + fak * S[1][7][3];
                S[5][7][3] = gma1 * S[2][7][3] + fak * (S[2][4][3] + S[2][7][1]);
                S[6][7][3] = gma2 * S[2][7][3] + fak * S[2][3][3];
                S[8][7][4] = gma0 * S[2][7][4] + fak * S[1][7][4];
                S[5][7][4] = gma1 * S[2][7][4] + fak * S[2][4][4];
                S[6][7][4] = gma2 * S[2][7][4] + fak * (S[2][3][4] + S[2][7][1]);
                S[8][8][2] = gma0 * S[2][8][2] + fak * (S[1][8][2] + S[2][2][2] + S[2][2][2] + S[2][8][1]);
                S[5][8][2] = gma1 * S[2][8][2];
                S[6][8][2] = gma2 * S[2][8][2];
                S[8][8][3] = gma0 * S[2][8][3] + fak * (S[1][8][3] + S[2][2][3] + S[2][2][3]);
                S[5][8][3] = gma1 * S[2][8][3] + fak * S[2][8][1];
                S[6][8][3] = gma2 * S[2][8][3];
                S[8][8][4] = gma0 * S[2][8][4] + fak * (S[1][8][4] + S[2][2][4] + S[2][2][4]);
                S[5][8][4] = gma1 * S[2][8][4];
                S[6][8][4] = gma2 * S[2][8][4] + fak * S[2][8][1];
                S[8][9][2] = gma0 * S[2][9][2] + fak * (S[1][9][2] + S[2][9][1]);
                S[5][9][2] = gma1 * S[2][9][2] + fak2 * S[2][3][2];
                S[6][9][2] = gma2 * S[2][9][2];
                S[8][9][3] = gma0 * S[2][9][3] + fak * S[1][9][3];
                S[5][9][3] = gma1 * S[2][9][3] + fak * (S[2][3][3] + S[2][3][3] + S[2][9][1]);
                S[6][9][3] = gma2 * S[2][9][3];
                S[8][9][4] = gma0 * S[2][9][4] + fak * S[1][9][4];
                S[5][9][4] = gma1 * S[2][9][4] + fak * (2 * S[2][3][4]);
                S[6][9][4] = gma2 * S[2][9][4] + fak * S[2][9][1];
                S[8][10][2] = gma0 * S[2][10][2] + fak * (S[1][10][2] + S[2][10][1]);
                S[5][10][2] = gma1 * S[2][10][2];
                S[6][10][2] = gma2 * S[2][10][2] + fak * (2 * S[2][4][2]);
                S[8][10][3] = gma0 * S[2][10][3] + fak * S[1][10][3];
                S[5][10][3] = gma1 * S[2][10][3] + fak * S[2][10][1];
                S[6][10][3] = gma2 * S[2][10][3] + fak * (2 * S[2][4][3]);
                S[8][10][4] = gma0 * S[2][10][4] + fak * S[1][10][4];
                S[5][10][4] = gma1 * S[2][10][4];
                S[6][10][4] = gma2 * S[2][10][4] + fak * (S[2][4][4] + S[2][4][4] + S[2][10][1]);
                S[9][5][2] = gma1 * S[3][5][2] + fak * (S[1][5][2] + S[3][2][2]);
                S[7][5][2] = gma2 * S[3][5][2];
                S[9][5][3] = gma1 * S[3][5][3] + fak * (S[1][5][3] + S[3][2][3] + S[3][5][1]);
                S[7][5][3] = gma2 * S[3][5][3];
                S[9][5][4] = gma1 * S[3][5][4] + fak * (S[1][5][4] + S[3][2][4]);
                S[7][5][4] = gma2 * S[3][5][4] + fak * S[3][5][1];
                S[9][6][2] = gma1 * S[3][6][2] + fak * S[1][6][2];
                S[7][6][2] = gma2 * S[3][6][2] + fak * S[3][2][2];
                S[9][6][3] = gma1 * S[3][6][3] + fak * (S[1][6][3] + S[3][6][1]);
                S[7][6][3] = gma2 * S[3][6][3] + fak * S[3][2][3];
                S[9][6][4] = gma1 * S[3][6][4] + fak * S[1][6][4];
                S[7][6][4] = gma2 * S[3][6][4] + fak * (S[3][2][4] + S[3][6][1]);
                S[9][7][2] = gma1 * S[3][7][2] + fak * (S[1][7][2] + S[3][4][2]);
                S[7][7][2] = gma2 * S[3][7][2] + fak * S[3][3][2];
                S[9][7][3] = gma1 * S[3][7][3] + fak * (S[1][7][3] + S[3][4][3] + S[3][7][1]);
                S[7][7][3] = gma2 * S[3][7][3] + fak * S[3][3][3];
                S[9][7][4] = gma1 * S[3][7][4] + fak * (S[1][7][4] + S[3][4][4]);
                S[7][7][4] = gma2 * S[3][7][4] + fak * (S[3][3][4] + S[3][7][1]);
                S[9][8][2] = gma1 * S[3][8][2] + fak * S[1][8][2];
                S[7][8][2] = gma2 * S[3][8][2];
                S[9][8][3] = gma1 * S[3][8][3] + fak * (S[1][8][3] + S[3][8][1]);
                S[7][8][3] = gma2 * S[3][8][3];
                S[9][8][4] = gma1 * S[3][8][4] + fak * S[1][8][4];
                S[7][8][4] = gma2 * S[3][8][4] + fak * S[3][8][1];
                S[9][9][2] = gma1 * S[3][9][2] + fak * (S[1][9][2] + S[3][3][2] + S[3][3][2]);
                S[7][9][2] = gma2 * S[3][9][2];
                S[9][9][3] = gma1 * S[3][9][3] + fak * (S[1][9][3] + S[3][3][3] + S[3][3][3] + S[3][9][1]);
                S[7][9][3] = gma2 * S[3][9][3];
                S[9][9][4] = gma1 * S[3][9][4] + fak * (S[1][9][4] + S[3][3][4] + S[3][3][4]);
                S[7][9][4] = gma2 * S[3][9][4] + fak * S[3][9][1];
                S[9][10][2] = gma1 * S[3][10][2] + fak * S[1][10][2];
                S[7][10][2] = gma2 * S[3][10][2] + fak2 * S[3][4][2];
                S[9][10][3] = gma1 * S[3][10][3] + fak * (S[1][10][3] + S[3][10][1]);
                S[7][10][3] = gma2 * S[3][10][3] + fak2 * S[3][4][3];
                S[9][10][4] = gma1 * S[3][10][4] + fak * S[1][10][4];
                S[7][10][4] = gma2 * S[3][10][4] + fak * (S[3][4][4] + S[3][4][4] + S[3][10][1]);
                S[10][5][2] = gma2 * S[4][5][2] + fak * S[1][5][2];
                S[10][5][3] = gma2 * S[4][5][3] + fak * S[1][5][3];
                S[10][5][4] = gma2 * S[4][5][4] + fak * (S[1][5][4] + S[4][5][1]);
                S[10][6][2] = gma2 * S[4][6][2] + fak * (S[1][6][2] + S[4][2][2]);
                S[10][6][3] = gma2 * S[4][6][3] + fak * (S[1][6][3] + S[4][2][3]);
                S[10][6][4] = gma2 * S[4][6][4] + fak * (S[1][6][4] + S[4][2][4] + S[4][6][1]);
                S[10][7][2] = gma2 * S[4][7][2] + fak * (S[1][7][2] + S[4][3][2]);
                S[10][7][3] = gma2 * S[4][7][3] + fak * (S[1][7][3] + S[4][3][3]);
                S[10][7][4] = gma2 * S[4][7][4] + fak * (S[1][7][4] + S[4][3][4] + S[4][7][1]);
                S[10][8][2] = gma2 * S[4][8][2] + fak * S[1][8][2];
                S[10][8][3] = gma2 * S[4][8][3] + fak * S[1][8][3];
                S[10][8][4] = gma2 * S[4][8][4] + fak * (S[1][8][4] + S[4][8][1]);
                S[10][9][2] = gma2 * S[4][9][2] + fak * S[1][9][2];
                S[10][9][3] = gma2 * S[4][9][3] + fak * S[1][9][3];
                S[10][9][4] = gma2 * S[4][9][4] + fak * (S[1][9][4] + S[4][9][1]);
                S[10][10][2] = gma2 * S[4][10][2] + fak * (S[1][10][2] + S[4][4][2] + S[4][4][2]);
                S[10][10][3] = gma2 * S[4][10][3] + fak * (S[1][10][3] + S[4][4][3] + S[4][4][3]);
                S[10][10][4] = gma2 * S[4][10][4] + fak * (S[1][10][4] + S[4][4][4] + S[4][4][4] + S[4][10][1]);
            }

            if (_lmax_alpha >= 2 && _lmax_gw >= 0 && _lmax_gamma >= 2) {
                S[8][1][5] = gma0 * S[2][1][5] + fak * (S[1][1][5] + S[2][1][3]);
                S[5][1][5] = gma1 * S[2][1][5] + fak * S[2][1][2];
                S[6][1][5] = gma2 * S[2][1][5];
                S[8][1][6] = gma0 * S[2][1][6] + fak * (S[1][1][6] + S[2][1][4]);
                S[5][1][6] = gma1 * S[2][1][6];
                S[6][1][6] = gma2 * S[2][1][6] + fak * S[2][1][2];
                S[8][1][7] = gma0 * S[2][1][7] + fak * S[1][1][7];
                S[5][1][7] = gma1 * S[2][1][7] + fak * S[2][1][4];
                S[6][1][7] = gma2 * S[2][1][7] + fak * S[2][1][3];
                S[8][1][8] = gma0 * S[2][1][8] + fak * (S[1][1][8] + S[2][1][2] + S[2][1][2]);
                S[5][1][8] = gma1 * S[2][1][8];
                S[6][1][8] = gma2 * S[2][1][8];
                S[8][1][9] = gma0 * S[2][1][9] + fak * S[1][1][9];
                S[5][1][9] = gma1 * S[2][1][9] + fak2 * S[2][1][3];
                S[6][1][9] = gma2 * S[2][1][9];
                S[8][1][10] = gma0 * S[2][1][10] + fak * S[1][1][10];
                S[5][1][10] = gma1 * S[2][1][10];
                S[6][1][10] = gma2 * S[2][1][10] + fak2 * S[2][1][4];
                S[9][1][5] = gma1 * S[3][1][5] + fak * (S[1][1][5] + S[3][1][2]);
                S[7][1][5] = gma2 * S[3][1][5];
                S[9][1][6] = gma1 * S[3][1][6] + fak * S[1][1][6];
                S[7][1][6] = gma2 * S[3][1][6] + fak * S[3][1][2];
                S[9][1][7] = gma1 * S[3][1][7] + fak * (S[1][1][7] + S[3][1][4]);
                S[7][1][7] = gma2 * S[3][1][7] + fak * S[3][1][3];
                S[9][1][8] = gma1 * S[3][1][8] + fak * S[1][1][8];
                S[7][1][8] = gma2 * S[3][1][8];
                S[9][1][9] = gma1 * S[3][1][9] + fak * (S[1][1][9] + S[3][1][3] + S[3][1][3]);
                S[7][1][9] = gma2 * S[3][1][9];
                S[9][1][10] = gma1 * S[3][1][10] + fak * S[1][1][10];
                S[7][1][10] = gma2 * S[3][1][10] + fak2 * S[3][1][4];
                S[10][1][5] = gma2 * S[4][1][5] + fak * S[1][1][5];
                S[10][1][6] = gma2 * S[4][1][6] + fak * (S[1][1][6] + S[4][1][2]);
                S[10][1][7] = gma2 * S[4][1][7] + fak * (S[1][1][7] + S[4][1][3]);
                S[10][1][8] = gma2 * S[4][1][8] + fak * S[1][1][8];
                S[10][1][9] = gma2 * S[4][1][9] + fak * S[1][1][9];
                S[10][1][10] = gma2 * S[4][1][10] + fak * (S[1][1][10] + S[4][1][4] + S[4][1][4]);
            }

            if (_lmax_alpha >= 2 && _lmax_gw >= 1 && _lmax_gamma >= 2) {
                S[8][2][5] = gma0 * S[2][2][5] + fak * (S[1][2][5] + S[2][1][5] + S[2][2][3]);
                S[5][2][5] = gma1 * S[2][2][5] + fak * S[2][2][2];
                S[6][2][5] = gma2 * S[2][2][5];
                S[8][2][6] = gma0 * S[2][2][6] + fak * (S[1][2][6] + S[2][1][6] + S[2][2][4]);
                S[5][2][6] = gma1 * S[2][2][6];
                S[6][2][6] = gma2 * S[2][2][6] + fak * S[2][2][2];
                S[8][2][7] = gma0 * S[2][2][7] + fak * (S[1][2][7] + S[2][1][7]);
                S[5][2][7] = gma1 * S[2][2][7] + fak * S[2][2][4];
                S[6][2][7] = gma2 * S[2][2][7] + fak * S[2][2][3];
                S[8][2][8] = gma0 * S[2][2][8] + fak * (S[1][2][8] + S[2][1][8] + S[2][2][2] + S[2][2][2]);
                S[5][2][8] = gma1 * S[2][2][8];
                S[6][2][8] = gma2 * S[2][2][8];
                S[8][2][9] = gma0 * S[2][2][9] + fak * (S[1][2][9] + S[2][1][9]);
                S[5][2][9] = gma1 * S[2][2][9] + fak2 * S[2][2][3];
                S[6][2][9] = gma2 * S[2][2][9];
                S[8][2][10] = gma0 * S[2][2][10] + fak * (S[1][2][10] + S[2][1][10]);
                S[5][2][10] = gma1 * S[2][2][10];
                S[6][2][10] = gma2 * S[2][2][10] + fak2 * S[2][2][4];
                S[8][3][5] = gma0 * S[2][3][5] + fak * (S[1][3][5] + S[2][3][3]);
                S[5][3][5] = gma1 * S[2][3][5] + fak * (S[2][1][5] + S[2][3][2]);
                S[6][3][5] = gma2 * S[2][3][5];
                S[8][3][6] = gma0 * S[2][3][6] + fak * (S[1][3][6] + S[2][3][4]);
                S[5][3][6] = gma1 * S[2][3][6] + fak * S[2][1][6];
                S[6][3][6] = gma2 * S[2][3][6] + fak * S[2][3][2];
                S[8][3][7] = gma0 * S[2][3][7] + fak * S[1][3][7];
                S[5][3][7] = gma1 * S[2][3][7] + fak * (S[2][1][7] + S[2][3][4]);
                S[6][3][7] = gma2 * S[2][3][7] + fak * S[2][3][3];
                S[8][3][8] = gma0 * S[2][3][8] + fak * (S[1][3][8] + S[2][3][2] + S[2][3][2]);
                S[5][3][8] = gma1 * S[2][3][8] + fak * S[2][1][8];
                S[6][3][8] = gma2 * S[2][3][8];
                S[8][3][9] = gma0 * S[2][3][9] + fak * S[1][3][9];
                S[5][3][9] = gma1 * S[2][3][9] + fak * (S[2][1][9] + S[2][3][3] + S[2][3][3]);
                S[6][3][9] = gma2 * S[2][3][9];
                S[8][3][10] = gma0 * S[2][3][10] + fak * S[1][3][10];
                S[5][3][10] = gma1 * S[2][3][10] + fak * S[2][1][10];
                S[6][3][10] = gma2 * S[2][3][10] + fak2 * S[2][3][4];
                S[8][4][5] = gma0 * S[2][4][5] + fak * (S[1][4][5] + S[2][4][3]);
                S[5][4][5] = gma1 * S[2][4][5] + fak * S[2][4][2];
                S[6][4][5] = gma2 * S[2][4][5] + fak * S[2][1][5];
                S[8][4][6] = gma0 * S[2][4][6] + fak * (S[1][4][6] + S[2][4][4]);
                S[5][4][6] = gma1 * S[2][4][6];
                S[6][4][6] = gma2 * S[2][4][6] + fak * (S[2][1][6] + S[2][4][2]);
                S[8][4][7] = gma0 * S[2][4][7] + fak * S[1][4][7];
                S[5][4][7] = gma1 * S[2][4][7] + fak * S[2][4][4];
                S[6][4][7] = gma2 * S[2][4][7] + fak * (S[2][1][7] + S[2][4][3]);
                S[8][4][8] = gma0 * S[2][4][8] + fak * (S[1][4][8] + S[2][4][2] + S[2][4][2]);
                S[5][4][8] = gma1 * S[2][4][8];
                S[6][4][8] = gma2 * S[2][4][8] + fak * S[2][1][8];
                S[8][4][9] = gma0 * S[2][4][9] + fak * S[1][4][9];
                S[5][4][9] = gma1 * S[2][4][9] + fak2 * S[2][4][3];
                S[6][4][9] = gma2 * S[2][4][9] + fak * S[2][1][9];
                S[8][4][10] = gma0 * S[2][4][10] + fak * S[1][4][10];
                S[5][4][10] = gma1 * S[2][4][10];
                S[6][4][10] = gma2 * S[2][4][10] + fak * (S[2][1][10] + S[2][4][4] + S[2][4][4]);
                S[9][2][5] = gma1 * S[3][2][5] + fak * (S[1][2][5] + S[3][2][2]);
                S[7][2][5] = gma2 * S[3][2][5];
                S[9][2][6] = gma1 * S[3][2][6] + fak * S[1][2][6];
                S[7][2][6] = gma2 * S[3][2][6] + fak * S[3][2][2];
                S[9][2][7] = gma1 * S[3][2][7] + fak * (S[1][2][7] + S[3][2][4]);
                S[7][2][7] = gma2 * S[3][2][7] + fak * S[3][2][3];
                S[9][2][8] = gma1 * S[3][2][8] + fak * S[1][2][8];
                S[7][2][8] = gma2 * S[3][2][8];
                S[9][2][9] = gma1 * S[3][2][9] + fak * (S[1][2][9] + S[3][2][3] + S[3][2][3]);
                S[7][2][9] = gma2 * S[3][2][9];
                S[9][2][10] = gma1 * S[3][2][10] + fak * S[1][2][10];
                S[7][2][10] = gma2 * S[3][2][10] + fak2 * S[3][2][4];
                S[9][3][5] = gma1 * S[3][3][5] + fak * (S[1][3][5] + S[3][1][5] + S[3][3][2]);
                S[7][3][5] = gma2 * S[3][3][5];
                S[9][3][6] = gma1 * S[3][3][6] + fak * (S[1][3][6] + S[3][1][6]);
                S[7][3][6] = gma2 * S[3][3][6] + fak * S[3][3][2];
                S[9][3][7] = gma1 * S[3][3][7] + fak * (S[1][3][7] + S[3][1][7] + S[3][3][4]);
                S[7][3][7] = gma2 * S[3][3][7] + fak * S[3][3][3];
                S[9][3][8] = gma1 * S[3][3][8] + fak * (S[1][3][8] + S[3][1][8]);
                S[7][3][8] = gma2 * S[3][3][8];
                S[9][3][9] = gma1 * S[3][3][9] + fak * (S[1][3][9] + S[3][1][9] + S[3][3][3] + S[3][3][3]);
                S[7][3][9] = gma2 * S[3][3][9];
                S[9][3][10] = gma1 * S[3][3][10] + fak * (S[1][3][10] + S[3][1][10]);
                S[7][3][10] = gma2 * S[3][3][10] + fak2 * S[3][3][4];
                S[9][4][5] = gma1 * S[3][4][5] + fak * (S[1][4][5] + S[3][4][2]);
                S[7][4][5] = gma2 * S[3][4][5] + fak * S[3][1][5];
                S[9][4][6] = gma1 * S[3][4][6] + fak * S[1][4][6];
                S[7][4][6] = gma2 * S[3][4][6] + fak * (S[3][1][6] + S[3][4][2]);
                S[9][4][7] = gma1 * S[3][4][7] + fak * (S[1][4][7] + S[3][4][4]);
                S[7][4][7] = gma2 * S[3][4][7] + fak * (S[3][1][7] + S[3][4][3]);
                S[9][4][8] = gma1 * S[3][4][8] + fak * S[1][4][8];
                S[7][4][8] = gma2 * S[3][4][8] + fak * S[3][1][8];
                S[9][4][9] = gma1 * S[3][4][9] + fak * (S[1][4][9] + S[3][4][3] + S[3][4][3]);
                S[7][4][9] = gma2 * S[3][4][9] + fak * S[3][1][9];
                S[9][4][10] = gma1 * S[3][4][10] + fak * S[1][4][10];
                S[7][4][10] = gma2 * S[3][4][10] + fak * (S[3][1][10] + S[3][4][4] + S[3][4][4]);
                S[10][2][5] = gma2 * S[4][2][5] + fak * S[1][2][5];
                S[10][2][6] = gma2 * S[4][2][6] + fak * (S[1][2][6] + S[4][2][2]);
                S[10][2][7] = gma2 * S[4][2][7] + fak * (S[1][2][7] + S[4][2][3]);
                S[10][2][8] = gma2 * S[4][2][8] + fak * S[1][2][8];
                S[10][2][9] = gma2 * S[4][2][9] + fak * S[1][2][9];
                S[10][2][10] = gma2 * S[4][2][10] + fak * (S[1][2][10] + S[4][2][4] + S[4][2][4]);
                S[10][3][5] = gma2 * S[4][3][5] + fak * S[1][3][5];
                S[10][3][6] = gma2 * S[4][3][6] + fak * (S[1][3][6] + S[4][3][2]);
                S[10][3][7] = gma2 * S[4][3][7] + fak * (S[1][3][7] + S[4][3][3]);
                S[10][3][8] = gma2 * S[4][3][8] + fak * S[1][3][8];
                S[10][3][9] = gma2 * S[4][3][9] + fak * S[1][3][9];
                S[10][3][10] = gma2 * S[4][3][10] + fak * (S[1][3][10] + S[4][3][4] + S[4][3][4]);
                S[10][4][5] = gma2 * S[4][4][5] + fak * (S[1][4][5] + S[4][1][5]);
                S[10][4][6] = gma2 * S[4][4][6] + fak * (S[1][4][6] + S[4][1][6] + S[4][4][2]);
                S[10][4][7] = gma2 * S[4][4][7] + fak * (S[1][4][7] + S[4][1][7] + S[4][4][3]);
                S[10][4][8] = gma2 * S[4][4][8] + fak * (S[1][4][8] + S[4][1][8]);
                S[10][4][9] = gma2 * S[4][4][9] + fak * (S[1][4][9] + S[4][1][9]);
                S[10][4][10] = gma2 * S[4][4][10] + fak * (S[1][4][10] + S[4][1][10] + S[4][4][4] + S[4][4][4]);
            }

            if (_lmax_alpha >= 2 && _lmax_gw >= 2 && _lmax_gamma >= 2) {
                S[8][5][5] = gma0 * S[2][5][5] + fak * (S[1][5][5] + S[2][3][5] + S[2][5][3]);
                S[5][5][5] = gma1 * S[2][5][5] + fak * (S[2][2][5] + S[2][5][2]);
                S[6][5][5] = gma2 * S[2][5][5];
                S[8][5][6] = gma0 * S[2][5][6] + fak * (S[1][5][6] + S[2][3][6] + S[2][5][4]);
                S[5][5][6] = gma1 * S[2][5][6] + fak * S[2][2][6];
                S[6][5][6] = gma2 * S[2][5][6] + fak * S[2][5][2];
                S[8][5][7] = gma0 * S[2][5][7] + fak * (S[1][5][7] + S[2][3][7]);
                S[5][5][7] = gma1 * S[2][5][7] + fak * (S[2][2][7] + S[2][5][4]);
                S[6][5][7] = gma2 * S[2][5][7] + fak * S[2][5][3];
                S[8][5][8] = gma0 * S[2][5][8] + fak * (S[1][5][8] + S[2][3][8] + S[2][5][2] + S[2][5][2]);
                S[5][5][8] = gma1 * S[2][5][8] + fak * S[2][2][8];
                S[6][5][8] = gma2 * S[2][5][8];
                S[8][5][9] = gma0 * S[2][5][9] + fak * (S[1][5][9] + S[2][3][9]);
                S[5][5][9] = gma1 * S[2][5][9] + fak * (S[2][2][9] + S[2][5][3] + S[2][5][3]);
                S[6][5][9] = gma2 * S[2][5][9];
                S[8][5][10] = gma0 * S[2][5][10] + fak * (S[1][5][10] + S[2][3][10]);
                S[5][5][10] = gma1 * S[2][5][10] + fak * S[2][2][10];
                S[6][5][10] = gma2 * S[2][5][10] + fak2 * S[2][5][4];
                S[8][6][5] = gma0 * S[2][6][5] + fak * (S[1][6][5] + S[2][4][5] + S[2][6][3]);
                S[5][6][5] = gma1 * S[2][6][5] + fak * S[2][6][2];
                S[6][6][5] = gma2 * S[2][6][5] + fak * S[2][2][5];
                S[8][6][6] = gma0 * S[2][6][6] + fak * (S[1][6][6] + S[2][4][6] + S[2][6][4]);
                S[5][6][6] = gma1 * S[2][6][6];
                S[6][6][6] = gma2 * S[2][6][6] + fak * (S[2][2][6] + S[2][6][2]);
                S[8][6][7] = gma0 * S[2][6][7] + fak * (S[1][6][7] + S[2][4][7]);
                S[5][6][7] = gma1 * S[2][6][7] + fak * S[2][6][4];
                S[6][6][7] = gma2 * S[2][6][7] + fak * (S[2][2][7] + S[2][6][3]);
                S[8][6][8] = gma0 * S[2][6][8] + fak * (S[1][6][8] + S[2][4][8] + S[2][6][2] + S[2][6][2]);
                S[5][6][8] = gma1 * S[2][6][8];
                S[6][6][8] = gma2 * S[2][6][8] + fak * S[2][2][8];
                S[8][6][9] = gma0 * S[2][6][9] + fak * (S[1][6][9] + S[2][4][9]);
                S[5][6][9] = gma1 * S[2][6][9] + fak2 * S[2][6][3];
                S[6][6][9] = gma2 * S[2][6][9] + fak * S[2][2][9];
                S[8][6][10] = gma0 * S[2][6][10] + fak * (S[1][6][10] + S[2][4][10]);
                S[5][6][10] = gma1 * S[2][6][10];
                S[6][6][10] = gma2 * S[2][6][10] + fak * (S[2][2][10] + S[2][6][4] + S[2][6][4]);
                S[8][7][5] = gma0 * S[2][7][5] + fak * (S[1][7][5] + S[2][7][3]);
                S[5][7][5] = gma1 * S[2][7][5] + fak * (S[2][4][5] + S[2][7][2]);
                S[6][7][5] = gma2 * S[2][7][5] + fak * S[2][3][5];
                S[8][7][6] = gma0 * S[2][7][6] + fak * (S[1][7][6] + S[2][7][4]);
                S[5][7][6] = gma1 * S[2][7][6] + fak * S[2][4][6];
                S[6][7][6] = gma2 * S[2][7][6] + fak * (S[2][3][6] + S[2][7][2]);
                S[8][7][7] = gma0 * S[2][7][7] + fak * S[1][7][7];
                S[5][7][7] = gma1 * S[2][7][7] + fak * (S[2][4][7] + S[2][7][4]);
                S[6][7][7] = gma2 * S[2][7][7] + fak * (S[2][3][7] + S[2][7][3]);
                S[8][7][8] = gma0 * S[2][7][8] + fak * (S[1][7][8] + S[2][7][2] + S[2][7][2]);
                S[5][7][8] = gma1 * S[2][7][8] + fak * S[2][4][8];
                S[6][7][8] = gma2 * S[2][7][8] + fak * S[2][3][8];
                S[8][7][9] = gma0 * S[2][7][9] + fak * S[1][7][9];
                S[5][7][9] = gma1 * S[2][7][9] + fak * (S[2][4][9] + S[2][7][3] + S[2][7][3]);
                S[6][7][9] = gma2 * S[2][7][9] + fak * S[2][3][9];
                S[8][7][10] = gma0 * S[2][7][10] + fak * S[1][7][10];
                S[5][7][10] = gma1 * S[2][7][10] + fak * S[2][4][10];
                S[6][7][10] = gma2 * S[2][7][10] + fak * (S[2][3][10] + S[2][7][4] + S[2][7][4]);
                S[8][8][5] = gma0 * S[2][8][5] + fak * (S[1][8][5] + S[2][2][5] + S[2][2][5] + S[2][8][3]);
                S[5][8][5] = gma1 * S[2][8][5] + fak * S[2][8][2];
                S[6][8][5] = gma2 * S[2][8][5];
                S[8][8][6] = gma0 * S[2][8][6] + fak * (S[1][8][6] + S[2][2][6] + S[2][2][6] + S[2][8][4]);
                S[5][8][6] = gma1 * S[2][8][6];
                S[6][8][6] = gma2 * S[2][8][6] + fak * S[2][8][2];
                S[8][8][7] = gma0 * S[2][8][7] + fak * (S[1][8][7] + S[2][2][7] + S[2][2][7]);
                S[5][8][7] = gma1 * S[2][8][7] + fak * S[2][8][4];
                S[6][8][7] = gma2 * S[2][8][7] + fak * S[2][8][3];
                S[8][8][8] = gma0 * S[2][8][8] + fak * (S[1][8][8] + S[2][2][8] + S[2][2][8] + S[2][8][2] + S[2][8][2]);
                S[5][8][8] = gma1 * S[2][8][8];
                S[6][8][8] = gma2 * S[2][8][8];
                S[8][8][9] = gma0 * S[2][8][9] + fak * (S[1][8][9] + S[2][2][9] + S[2][2][9]);
                S[5][8][9] = gma1 * S[2][8][9] + fak2 * S[2][8][3];
                S[6][8][9] = gma2 * S[2][8][9];
                S[8][8][10] = gma0 * S[2][8][10] + fak * (S[1][8][10] + S[2][2][10] + S[2][2][10]);
                S[5][8][10] = gma1 * S[2][8][10];
                S[6][8][10] = gma2 * S[2][8][10] + fak2 * S[2][8][4];
                S[8][9][5] = gma0 * S[2][9][5] + fak * (S[1][9][5] + S[2][9][3]);
                S[5][9][5] = gma1 * S[2][9][5] + fak * (S[2][3][5] + S[2][3][5] + S[2][9][2]);
                S[6][9][5] = gma2 * S[2][9][5];
                S[8][9][6] = gma0 * S[2][9][6] + fak * (S[1][9][6] + S[2][9][4]);
                S[5][9][6] = gma1 * S[2][9][6] + fak2 * S[2][3][6];
                S[6][9][6] = gma2 * S[2][9][6] + fak * S[2][9][2];
                S[8][9][7] = gma0 * S[2][9][7] + fak * S[1][9][7];
                S[5][9][7] = gma1 * S[2][9][7] + fak * (S[2][3][7] + S[2][3][7] + S[2][9][4]);
                S[6][9][7] = gma2 * S[2][9][7] + fak * S[2][9][3];
                S[8][9][8] = gma0 * S[2][9][8] + fak * (S[1][9][8] + S[2][9][2] + S[2][9][2]);
                S[5][9][8] = gma1 * S[2][9][8] + fak2 * S[2][3][8];
                S[6][9][8] = gma2 * S[2][9][8];
                S[8][9][9] = gma0 * S[2][9][9] + fak * S[1][9][9];
                S[5][9][9] = gma1 * S[2][9][9] + fak2 * (S[2][3][9] + S[2][9][3]);
                S[6][9][9] = gma2 * S[2][9][9];
                S[8][9][10] = gma0 * S[2][9][10] + fak * S[1][9][10];
                S[5][9][10] = gma1 * S[2][9][10] + fak2 * S[2][3][10];
                S[6][9][10] = gma2 * S[2][9][10] + fak2 * S[2][9][4];
                S[8][10][5] = gma0 * S[2][10][5] + fak * (S[1][10][5] + S[2][10][3]);
                S[5][10][5] = gma1 * S[2][10][5] + fak * S[2][10][2];
                S[6][10][5] = gma2 * S[2][10][5] + fak2 * S[2][4][5];
                S[8][10][6] = gma0 * S[2][10][6] + fak * (S[1][10][6] + S[2][10][4]);
                S[5][10][6] = gma1 * S[2][10][6];
                S[6][10][6] = gma2 * S[2][10][6] + fak * (S[2][4][6] + S[2][4][6] + S[2][10][2]);
                S[8][10][7] = gma0 * S[2][10][7] + fak * S[1][10][7];
                S[5][10][7] = gma1 * S[2][10][7] + fak * S[2][10][4];
                S[6][10][7] = gma2 * S[2][10][7] + fak * (S[2][4][7] + S[2][4][7] + S[2][10][3]);
                S[8][10][8] = gma0 * S[2][10][8] + fak * (S[1][10][8] + S[2][10][2] + S[2][10][2]);
                S[5][10][8] = gma1 * S[2][10][8];
                S[6][10][8] = gma2 * S[2][10][8] + fak2 * S[2][4][8];
                S[8][10][9] = gma0 * S[2][10][9] + fak * S[1][10][9];
                S[5][10][9] = gma1 * S[2][10][9] + fak2 * S[2][10][3];
                S[6][10][9] = gma2 * S[2][10][9] + fak2 * S[2][4][9];
                S[8][10][10] = gma0 * S[2][10][10] + fak * S[1][10][10];
                S[5][10][10] = gma1 * S[2][10][10];
                S[6][10][10] = gma2 * S[2][10][10] + fak2 * (S[2][4][10] + S[2][10][4]);
                S[9][5][5] = gma1 * S[3][5][5] + fak * (S[1][5][5] + S[3][2][5] + S[3][5][2]);
                S[7][5][5] = gma2 * S[3][5][5];
                S[9][5][6] = gma1 * S[3][5][6] + fak * (S[1][5][6] + S[3][2][6]);
                S[7][5][6] = gma2 * S[3][5][6] + fak * S[3][5][2];
                S[9][5][7] = gma1 * S[3][5][7] + fak * (S[1][5][7] + S[3][2][7] + S[3][5][4]);
                S[7][5][7] = gma2 * S[3][5][7] + fak * S[3][5][3];
                S[9][5][8] = gma1 * S[3][5][8] + fak * (S[1][5][8] + S[3][2][8]);
                S[7][5][8] = gma2 * S[3][5][8];
                S[9][5][9] = gma1 * S[3][5][9] + fak * (S[1][5][9] + S[3][2][9] + S[3][5][3] + S[3][5][3]);
                S[7][5][9] = gma2 * S[3][5][9];
                S[9][5][10] = gma1 * S[3][5][10] + fak * (S[1][5][10] + S[3][2][10]);
                S[7][5][10] = gma2 * S[3][5][10] + fak2 * S[3][5][4];
                S[9][6][5] = gma1 * S[3][6][5] + fak * (S[1][6][5] + S[3][6][2]);
                S[7][6][5] = gma2 * S[3][6][5] + fak * S[3][2][5];
                S[9][6][6] = gma1 * S[3][6][6] + fak * S[1][6][6];
                S[7][6][6] = gma2 * S[3][6][6] + fak * (S[3][2][6] + S[3][6][2]);
                S[9][6][7] = gma1 * S[3][6][7] + fak * (S[1][6][7] + S[3][6][4]);
                S[7][6][7] = gma2 * S[3][6][7] + fak * (S[3][2][7] + S[3][6][3]);
                S[9][6][8] = gma1 * S[3][6][8] + fak * S[1][6][8];
                S[7][6][8] = gma2 * S[3][6][8] + fak * S[3][2][8];
                S[9][6][9] = gma1 * S[3][6][9] + fak * (S[1][6][9] + S[3][6][3] + S[3][6][3]);
                S[7][6][9] = gma2 * S[3][6][9] + fak * S[3][2][9];
                S[9][6][10] = gma1 * S[3][6][10] + fak * S[1][6][10];
                S[7][6][10] = gma2 * S[3][6][10] + fak * (S[3][2][10] + S[3][6][4] + S[3][6][4]);
                S[9][7][5] = gma1 * S[3][7][5] + fak * (S[1][7][5] + S[3][4][5] + S[3][7][2]);
                S[7][7][5] = gma2 * S[3][7][5] + fak * S[3][3][5];
                S[9][7][6] = gma1 * S[3][7][6] + fak * (S[1][7][6] + S[3][4][6]);
                S[7][7][6] = gma2 * S[3][7][6] + fak * (S[3][3][6] + S[3][7][2]);
                S[9][7][7] = gma1 * S[3][7][7] + fak * (S[1][7][7] + S[3][4][7] + S[3][7][4]);
                S[7][7][7] = gma2 * S[3][7][7] + fak * (S[3][3][7] + S[3][7][3]);
                S[9][7][8] = gma1 * S[3][7][8] + fak * (S[1][7][8] + S[3][4][8]);
                S[7][7][8] = gma2 * S[3][7][8] + fak * S[3][3][8];
                S[9][7][9] = gma1 * S[3][7][9] + fak * (S[1][7][9] + S[3][4][9] + S[3][7][3] + S[3][7][3]);
                S[7][7][9] = gma2 * S[3][7][9] + fak * S[3][3][9];
                S[9][7][10] = gma1 * S[3][7][10] + fak * (S[1][7][10] + S[3][4][10]);
                S[7][7][10] = gma2 * S[3][7][10] + fak * (S[3][3][10] + S[3][7][4] + S[3][7][4]);
                S[9][8][5] = gma1 * S[3][8][5] + fak * (S[1][8][5] + S[3][8][2]);
                S[7][8][5] = gma2 * S[3][8][5];
                S[9][8][6] = gma1 * S[3][8][6] + fak * S[1][8][6];
                S[7][8][6] = gma2 * S[3][8][6] + fak * S[3][8][2];
                S[9][8][7] = gma1 * S[3][8][7] + fak * (S[1][8][7] + S[3][8][4]);
                S[7][8][7] = gma2 * S[3][8][7] + fak * S[3][8][3];
                S[9][8][8] = gma1 * S[3][8][8] + fak * S[1][8][8];
                S[7][8][8] = gma2 * S[3][8][8];
                S[9][8][9] = gma1 * S[3][8][9] + fak * (S[1][8][9] + S[3][8][3] + S[3][8][3]);
                S[7][8][9] = gma2 * S[3][8][9];
                S[9][8][10] = gma1 * S[3][8][10] + fak * S[1][8][10];
                S[7][8][10] = gma2 * S[3][8][10] + fak2 * S[3][8][4];
                S[9][9][5] = gma1 * S[3][9][5] + fak * (S[1][9][5] + S[3][3][5] + S[3][3][5] + S[3][9][2]);
                S[7][9][5] = gma2 * S[3][9][5];
                S[9][9][6] = gma1 * S[3][9][6] + fak * (S[1][9][6] + S[3][3][6] + S[3][3][6]);
                S[7][9][6] = gma2 * S[3][9][6] + fak * S[3][9][2];
                S[9][9][7] = gma1 * S[3][9][7] + fak * (S[1][9][7] + S[3][3][7] + S[3][3][7] + S[3][9][4]);
                S[7][9][7] = gma2 * S[3][9][7] + fak * S[3][9][3];
                S[9][9][8] = gma1 * S[3][9][8] + fak * (S[1][9][8] + S[3][3][8] + S[3][3][8]);
                S[7][9][8] = gma2 * S[3][9][8];
                S[9][9][9] = gma1 * S[3][9][9] + fak * (S[1][9][9] + S[3][3][9] + S[3][3][9] + S[3][9][3] + S[3][9][3]);
                S[7][9][9] = gma2 * S[3][9][9];
                S[9][9][10] = gma1 * S[3][9][10] + fak * (S[1][9][10] + S[3][3][10] + S[3][3][10]);
                S[7][9][10] = gma2 * S[3][9][10] + fak2 * S[3][9][4];
                S[9][10][5] = gma1 * S[3][10][5] + fak * (S[1][10][5] + S[3][10][2]);
                S[7][10][5] = gma2 * S[3][10][5] + fak2 * S[3][4][5];
                S[9][10][6] = gma1 * S[3][10][6] + fak * S[1][10][6];
                S[7][10][6] = gma2 * S[3][10][6] + fak * (S[3][4][6] + S[3][4][6] + S[3][10][2]);
                S[9][10][7] = gma1 * S[3][10][7] + fak * (S[1][10][7] + S[3][10][4]);
                S[7][10][7] = gma2 * S[3][10][7] + fak * (S[3][4][7] + S[3][4][7] + S[3][10][3]);
                S[9][10][8] = gma1 * S[3][10][8] + fak * S[1][10][8];
                S[7][10][8] = gma2 * S[3][10][8] + fak2 * S[3][4][8];
                S[9][10][9] = gma1 * S[3][10][9] + fak * (S[1][10][9] + S[3][10][3] + S[3][10][3]);
                S[7][10][9] = gma2 * S[3][10][9] + fak2 * S[3][4][9];
                S[9][10][10] = gma1 * S[3][10][10] + fak * S[1][10][10];
                S[7][10][10] = gma2 * S[3][10][10] + fak2 * (S[3][4][10] + S[3][10][4]);
                S[10][5][5] = gma2 * S[4][5][5] + fak * S[1][5][5];
                S[10][5][6] = gma2 * S[4][5][6] + fak * (S[1][5][6] + S[4][5][2]);
                S[10][5][7] = gma2 * S[4][5][7] + fak * (S[1][5][7] + S[4][5][3]);
                S[10][5][8] = gma2 * S[4][5][8] + fak * S[1][5][8];
                S[10][5][9] = gma2 * S[4][5][9] + fak * S[1][5][9];
                S[10][5][10] = gma2 * S[4][5][10] + fak * (S[1][5][10] + S[4][5][4] + S[4][5][4]);
                S[10][6][5] = gma2 * S[4][6][5] + fak * (S[1][6][5] + S[4][2][5]);
                S[10][6][6] = gma2 * S[4][6][6] + fak * (S[1][6][6] + S[4][2][6] + S[4][6][2]);
                S[10][6][7] = gma2 * S[4][6][7] + fak * (S[1][6][7] + S[4][2][7] + S[4][6][3]);
                S[10][6][8] = gma2 * S[4][6][8] + fak * (S[1][6][8] + S[4][2][8]);
                S[10][6][9] = gma2 * S[4][6][9] + fak * (S[1][6][9] + S[4][2][9]);
                S[10][6][10] = gma2 * S[4][6][10] + fak * (S[1][6][10] + S[4][2][10] + S[4][6][4] + S[4][6][4]);
                S[10][7][5] = gma2 * S[4][7][5] + fak * (S[1][7][5] + S[4][3][5]);
                S[10][7][6] = gma2 * S[4][7][6] + fak * (S[1][7][6] + S[4][3][6] + S[4][7][2]);
                S[10][7][7] = gma2 * S[4][7][7] + fak * (S[1][7][7] + S[4][3][7] + S[4][7][3]);
                S[10][7][8] = gma2 * S[4][7][8] + fak * (S[1][7][8] + S[4][3][8]);
                S[10][7][9] = gma2 * S[4][7][9] + fak * (S[1][7][9] + S[4][3][9]);
                S[10][7][10] = gma2 * S[4][7][10] + fak * (S[1][7][10] + S[4][3][10] + S[4][7][4] + S[4][7][4]);
                S[10][8][5] = gma2 * S[4][8][5] + fak * (1 * S[1][8][5]);
                S[10][8][6] = gma2 * S[4][8][6] + fak * (S[1][8][6] + S[4][8][2]);
                S[10][8][7] = gma2 * S[4][8][7] + fak * (S[1][8][7] + S[4][8][3]);
                S[10][8][8] = gma2 * S[4][8][8] + fak * S[1][8][8];
                S[10][8][9] = gma2 * S[4][8][9] + fak * S[1][8][9];
                S[10][8][10] = gma2 * S[4][8][10] + fak * (S[1][8][10] + S[4][8][4] + S[4][8][4]);
                S[10][9][5] = gma2 * S[4][9][5] + fak * S[1][9][5];
                S[10][9][6] = gma2 * S[4][9][6] + fak * (S[1][9][6] + S[4][9][2]);
                S[10][9][7] = gma2 * S[4][9][7] + fak * (S[1][9][7] + S[4][9][3]);
                S[10][9][8] = gma2 * S[4][9][8] + fak * S[1][9][8];
                S[10][9][9] = gma2 * S[4][9][9] + fak * S[1][9][9];
                S[10][9][10] = gma2 * S[4][9][10] + fak * (S[1][9][10] + S[4][9][4] + S[4][9][4]);
                S[10][10][5] = gma2 * S[4][10][5] + fak * (S[1][10][5] + S[4][4][5] + S[4][4][5]);
                S[10][10][6] = gma2 * S[4][10][6] + fak * (S[1][10][6] + S[4][4][6] + S[4][4][6] + S[4][10][2]);
                S[10][10][7] = gma2 * S[4][10][7] + fak * (S[1][10][7] + S[4][4][7] + S[4][4][7] + S[4][10][3]);
                S[10][10][8] = gma2 * S[4][10][8] + fak * (S[1][10][8] + S[4][4][8] + S[4][4][8]);
                S[10][10][9] = gma2 * S[4][10][9] + fak * (S[1][10][9] + S[4][4][9] + S[4][4][9]);
                S[10][10][10] = gma2 * S[4][10][10] + fak * (S[1][10][10] + S[4][4][10] + S[4][4][10] + S[4][10][4] + S[4][10][4]);
            }
            /*
             if ( _lmax_alpha >= 0 && _lmax_gw >= 3 && _lmax_gamma >= 0) {
                S[1][14][1] = gmc0*S[1][5][1] + fak * S[1][3][1]
                S[1][15][1] = gmc1*S[1][5][1] + fak * S[1][2][1]
                S[1][20][1] = gmc2*S[1][5][1]
                S[1][16][1] = gmc0*S[1][6][1] + fak * S[1][4][1]
                S[1][17][1] = gmc2*S[1][6][1] + fak * S[1][2][1]
                S[1][18][1] = gmc1*S[1][7][1] + fak * S[1][4][1]
                S[1][19][1] = gmc2*S[1][7][1] + fak * S[1][3][1]
                S[1][11][1] = gmc0*S[1][8][1] + fak2* S[1][2][1]
                S[1][12][1] = gmc1*S[1][9][1] + fak2* S[1][3][1]
                S[1][13][1] = gmc2*S[1][10][1] + fak2* S[1][4][1]
             }

             if ( _lmax_alpha >= 0 && _lmax_gw >= 3 && _lmax_gamma >= 1) {
                S[1][14][2] = gmc0*S[1][5][2] + fak * (S[1][3][2] +S[1][5][1] )
                S[1][15][2] = gmc1*S[1][5][2] + fak * S[1][2][2]
                S[1][20][2] = gmc2*S[1][5][2]
                S[1][14][3] = gmc0*S[1][5][3] + fak * S[1][3][3]
                S[1][15][3] = gmc1*S[1][5][3] + fak * (S[1][2][3] +S[1][5][1] )
                S[1][20][3] = gmc2*S[1][5][3]
                S[1][14][4] = gmc0*S[1][5][4] + fak * S[1][3][4]
                S[1][15][4] = gmc1*S[1][5][4] + fak * S[1][2][4]
                S[1][20][4] = gmc2*S[1][5][4] + fak * S[1][5][1]
                S[1][16][2] = gmc0*S[1][6][2] + fak * (S[1][4][2] +S[1][6][1] )
                S[1][17][2] = gmc2*S[1][6][2] + fak * S[1][2][2]
                S[1][16][3] = gmc0*S[1][6][3] + fak * S[1][4][3]
                S[1][17][3] = gmc2*S[1][6][3] + fak * S[1][2][3]
                S[1][16][4] = gmc0*S[1][6][4] + fak * S[1][4][4]
                S[1][17][4] = gmc2*S[1][6][4] + fak * (S[1][2][4] +S[1][6][1] )
                S[1][18][2] = gmc1*S[1][7][2] + fak * S[1][4][2] 
                S[1][19][2] = gmc2*S[1][7][2] + fak * S[1][3][2]
                S[1][18][3] = gmc1*S[1][7][3] + fak * (S[1][4][3] +S[1][7][1] )
                S[1][19][3] = gmc2*S[1][7][3] + fak * S[1][3][3]
                S[1][18][4] = gmc1*S[1][7][4] + fak * S[1][4][4]
                S[1][19][4] = gmc2*S[1][7][4] + fak * (S[1][3][4] +S[1][7][1] )
                S[1][11][2] = gmc0*S[1][8][2] + fak * (2.d0*S[1][2][2] +S[1][8][1] )
                S[1][11][3] = gmc0*S[1][8][3] + fak2* S[1][2][3]
                S[1][11][4] = gmc0*S[1][8][4] + fak2* S[1][2][4]
                S[1][12][2] = gmc1*S[1][9][2] + fak2* S[1][3][2]
                S[1][12][3] = gmc1*S[1][9][3] + fak * (2.d0*S[1][3][3] +S[1][9][1] )
                S[1][12][4] = gmc1*S[1][9][4] + fak2* S[1][3][4]
                S[1][13][2] = gmc2*S[1][10][2] + fak2* S[1][4][2]
                S[1][13][3] = gmc2*S[1][10][3] + fak2* S[1][4][3]
                S[1][13][4] = gmc2*S[1][10][4] + fak * (2.d0*S[1][4][4] +S[1][10][1] )
             }

             if ( _lmax_alpha >= 1 && _lmax_gw >= 3 && _lmax_gamma >= 0) {
                S[2][11][1] = gma0*S[1][11][1] + fak3* S[1][8][1]
                S[3][11][1] = gma1*S[1][11][1]
                S[4][11][1] = gma2*S[1][11][1]
                S[2][12][1] = gma0*S[1][12][1]
                S[3][12][1] = gma1*S[1][12][1] + fak3* S[1][9][1]
                S[4][12][1] = gma2*S[1][12][1]
                S[2][13][1] = gma0*S[1][13][1]
                S[3][13][1] = gma1*S[1][13][1]
                S[4][13][1] = gma2*S[1][13][1] + fak3* S[1][10][1]
                S[2][14][1] = gma0*S[1][14][1] + fak2* S[1][5][1]
                S[3][14][1] = gma1*S[1][14][1] + fak * S[1][8][1]
                S[4][14][1] = gma2*S[1][14][1]
                S[2][15][1] = gma0*S[1][15][1] + fak * S[1][9][1]
                S[3][15][1] = gma1*S[1][15][1] + fak2* S[1][5][1]
                S[4][15][1] = gma2*S[1][15][1]
                S[2][16][1] = gma0*S[1][16][1] + fak2* S[1][6][1]
                S[3][16][1] = gma1*S[1][16][1]
                S[4][16][1] = gma2*S[1][16][1] + fak * S[1][8][1]
                S[2][17][1] = gma0*S[1][17][1] + fak * S[1][10][1]
                S[3][17][1] = gma1*S[1][17][1]
                S[4][17][1] = gma2*S[1][17][1] + fak2* S[1][6][1]
                S[2][18][1] = gma0*S[1][18][1]
                S[3][18][1] = gma1*S[1][18][1] + fak2* S[1][7][1]
                S[4][18][1] = gma2*S[1][18][1] + fak * S[1][9][1]
                S[2][19][1] = gma0*S[1][19][1]
                S[3][19][1] = gma1*S[1][19][1] + fak * S[1][10][1]
                S[4][19][1] = gma2*S[1][19][1] + fak2* S[1][7][1]
                S[2][20][1] = gma0*S[1][20][1] + fak * S[1][7][1]
                S[3][20][1] = gma1*S[1][20][1] + fak * S[1][6][1]
                S[4][20][1] = gma2*S[1][20][1] + fak * S[1][5][1]
             }

             if ( _lmax_alpha >= 0 && _lmax_gw >= 3 && _lmax_gamma >= 2) {
                S[1][14][5] = gmc0*S[1][5][5] + fak * (S[1][3][5] +S[1][5][3] )
                S[1][15][5] = gmc1*S[1][5][5] + fak * (S[1][2][5] +S[1][5][2] )
                S[1][20][5] = gmc2*S[1][5][5]
                S[1][14][6] = gmc0*S[1][5][6] + fak * (S[1][3][6] +S[1][5][4] )
                S[1][15][6] = gmc1*S[1][5][6] + fak * S[1][2][6]
                S[1][20][6] = gmc2*S[1][5][6] + fak * S[1][5][2]
                S[1][14][7] = gmc0*S[1][5][7] + fak * S[1][3][7]
                S[1][15][7] = gmc1*S[1][5][7] + fak * (S[1][2][7] +S[1][5][4] )
                S[1][20][7] = gmc2*S[1][5][7] + fak * S[1][5][3]
                S[1][14][8] = gmc0*S[1][5][8] + fak * (S[1][3][8] +2.d0*S[1][5][2] )
                S[1][15][8] = gmc1*S[1][5][8] + fak * S[1][2][8]
                S[1][20][8] = gmc2*S[1][5][8]
                S[1][14][9] = gmc0*S[1][5][9] + fak * S[1][3][9]
                S[1][15][9] = gmc1*S[1][5][9] + fak * (S[1][2][9] +2.d0*S[1][5][3] )
                S[1][20][9] = gmc2*S[1][5][9]
                S[1][14][10] = gmc0*S[1][5][10] + fak * S[1][3][10]
                S[1][15][10] = gmc1*S[1][5][10] + fak * S[1][2][10]
                S[1][20][10] = gmc2*S[1][5][10] + fak2* S[1][5][4]
                S[1][16][5] = gmc0*S[1][6][5] + fak * (S[1][4][5] +S[1][6][3] )
                S[1][17][5] = gmc2*S[1][6][5] + fak * S[1][2][5]
                S[1][16][6] = gmc0*S[1][6][6] + fak * (S[1][4][6] +S[1][6][4] )
                S[1][17][6] = gmc2*S[1][6][6] + fak * (S[1][2][6] +S[1][6][2] )
                S[1][16][7] = gmc0*S[1][6][7] + fak * S[1][4][7]
                S[1][17][7] = gmc2*S[1][6][7] + fak * (S[1][2][7] +S[1][6][3] )
                S[1][16][8] = gmc0*S[1][6][8] + fak * (S[1][4][8] +2.d0*S[1][6][2] )
                S[1][17][8] = gmc2*S[1][6][8] + fak * S[1][2][8]
                S[1][16][9] = gmc0*S[1][6][9] + fak * S[1][4][9]
                S[1][17][9] = gmc2*S[1][6][9] + fak * S[1][2][9]
                S[1][16][10] = gmc0*S[1][6][10] + fak * S[1][4][10]
                S[1][17][10] = gmc2*S[1][6][10] + fak * (S[1][2][10] +2.d0*S[1][6][4] )
                S[1][18][5] = gmc1*S[1][7][5] + fak * (S[1][4][5] +S[1][7][2] )
                S[1][19][5] = gmc2*S[1][7][5] + fak * S[1][3][5]
                S[1][18][6] = gmc1*S[1][7][6] + fak * S[1][4][6]
                S[1][19][6] = gmc2*S[1][7][6] + fak * (S[1][3][6] +S[1][7][2] )
                S[1][18][7] = gmc1*S[1][7][7] + fak * (S[1][4][7] +S[1][7][4] )
                S[1][19][7] = gmc2*S[1][7][7] + fak * (S[1][3][7] +S[1][7][3] )
                S[1][18][8] = gmc1*S[1][7][8] + fak * S[1][4][8]
                S[1][19][8] = gmc2*S[1][7][8] + fak * S[1][3][8]
                S[1][18][9] = gmc1*S[1][7][9] + fak * (S[1][4][9] +2.d0*S[1][7][3] )
                S[1][19][9] = gmc2*S[1][7][9] + fak * S[1][3][9]
                S[1][18][10] = gmc1*S[1][7][10] + fak * S[1][4][10]
                S[1][19][10] = gmc2*S[1][7][10] + fak * (S[1][3][10] +2.d0*S[1][7][4] )
                S[1][11][5] = gmc0*S[1][8][5] + fak * (2.d0*S[1][2][5] +S[1][8][3] )
                S[1][11][6] = gmc0*S[1][8][6] + fak * (2.d0*S[1][2][6] +S[1][8][4] )
                S[1][11][7] = gmc0*S[1][8][7] + fak2* S[1][2][7]
                S[1][11][8] = gmc0*S[1][8][8] + fak2* (S[1][2][8] +S[1][8][2] )
                S[1][11][9] = gmc0*S[1][8][9] + fak2* S[1][2][9]
                S[1][11][10] = gmc0*S[1][8][10] + fak2* S[1][2][10]
                S[1][12][5] = gmc1*S[1][9][5] + fak * (2.d0*S[1][3][5] +S[1][9][2] )
                S[1][12][6] = gmc1*S[1][9][6] + fak2* S[1][3][6]
                S[1][12][7] = gmc1*S[1][9][7] + fak * (2.d0*S[1][3][7] +S[1][9][4] )
                S[1][12][8] = gmc1*S[1][9][8] + fak2* S[1][3][8]
                S[1][12][9] = gmc1*S[1][9][9] + fak2* (S[1][3][9] +S[1][9][3] )
                S[1][12][10] = gmc1*S[1][9][10] + fak2* S[1][3][10]
                S[1][13][5] = gmc2*S[1][10][5] + fak2* S[1][4][5]
                S[1][13][6] = gmc2*S[1][10][6] + fak * (2.d0*S[1][4][6] +S[1][10][2] )
                S[1][13][7] = gmc2*S[1][10][7] + fak * (2.d0*S[1][4][7] +S[1][10][3] )
                S[1][13][8] = gmc2*S[1][10][8] + fak2* S[1][4][8]
                S[1][13][9] = gmc2*S[1][10][9] + fak2* S[1][4][9]
                S[1][13][10] = gmc2*S[1][10][10] + fak2* (S[1][4][10] +S[1][10][4] )
             }

             if ( _lmax_alpha >= 1 && _lmax_gw >= 3 && _lmax_gamma >= 1) {
                S[2][11][2] = gma0*S[1][11][2] + fak * (3.d0*S[1][8][2] +S[1][11][1] )
                S[3][11][2] = gma1*S[1][11][2]
                S[4][11][2] = gma2*S[1][11][2]
                S[2][11][3] = gma0*S[1][11][3] + fak3* S[1][8][3]
                S[3][11][3] = gma1*S[1][11][3] + fak * S[1][11][1]
                S[4][11][3] = gma2*S[1][11][3]
                S[2][11][4] = gma0*S[1][11][4] + fak3* S[1][8][4]
                S[3][11][4] = gma1*S[1][11][4]
                S[4][11][4] = gma2*S[1][11][4] + fak * S[1][11][1]
                S[2][12][2] = gma0*S[1][12][2] + fak * S[1][12][1]
                S[3][12][2] = gma1*S[1][12][2] + fak3* S[1][9][2]
                S[4][12][2] = gma2*S[1][12][2]
                S[2][12][3] = gma0*S[1][12][3]
                S[3][12][3] = gma1*S[1][12][3] + fak * (3.d0*S[1][9][3] +S[1][12][1] )
                S[4][12][3] = gma2*S[1][12][3]
                S[2][12][4] = gma0*S[1][12][4]
                S[3][12][4] = gma1*S[1][12][4] + fak3* S[1][9][4]
                S[4][12][4] = gma2*S[1][12][4] + fak * S[1][12][1]
                S[2][13][2] = gma0*S[1][13][2] + fak * S[1][13][1]
                S[3][13][2] = gma1*S[1][13][2]
                S[4][13][2] = gma2*S[1][13][2] + fak3* S[1][10][2]
                S[2][13][3] = gma0*S[1][13][3]
                S[3][13][3] = gma1*S[1][13][3] + fak * S[1][13][1]
                S[4][13][3] = gma2*S[1][13][3] + fak3* S[1][10][3]
                S[2][13][4] = gma0*S[1][13][4]
                S[3][13][4] = gma1*S[1][13][4]
                S[4][13][4] = gma2*S[1][13][4] + fak * (3.d0*S[1][10][4] +S[1][13][1] )
                S[2][14][2] = gma0*S[1][14][2] + fak * (2.d0*S[1][5][2] +S[1][14][1] )
                S[3][14][2] = gma1*S[1][14][2] + fak * S[1][8][2]
                S[4][14][2] = gma2*S[1][14][2]
                S[2][14][3] = gma0*S[1][14][3] + fak2* S[1][5][3]
                S[3][14][3] = gma1*S[1][14][3] + fak * (S[1][8][3] +S[1][14][1] )
                S[4][14][3] = gma2*S[1][14][3]
                S[2][14][4] = gma0*S[1][14][4] + fak2* S[1][5][4]
                S[3][14][4] = gma1*S[1][14][4] + fak * S[1][8][4]
                S[4][14][4] = gma2*S[1][14][4] + fak * S[1][14][1]
                S[2][15][2] = gma0*S[1][15][2] + fak * (S[1][9][2] +S[1][15][1] )
                S[3][15][2] = gma1*S[1][15][2] + fak2* S[1][5][2]
                S[4][15][2] = gma2*S[1][15][2]
                S[2][15][3] = gma0*S[1][15][3] + fak * S[1][9][3]
                S[3][15][3] = gma1*S[1][15][3] + fak * (2.d0*S[1][5][3] +S[1][15][1] )
                S[4][15][3] = gma2*S[1][15][3]
                S[2][15][4] = gma0*S[1][15][4] + fak * S[1][9][4]
                S[3][15][4] = gma1*S[1][15][4] + fak2* S[1][5][4]
                S[4][15][4] = gma2*S[1][15][4] + fak * S[1][15][1]
                S[2][16][2] = gma0*S[1][16][2] + fak * (2.d0*S[1][6][2] +S[1][16][1] )
                S[3][16][2] = gma1*S[1][16][2]
                S[4][16][2] = gma2*S[1][16][2] + fak * S[1][8][2]
                S[2][16][3] = gma0*S[1][16][3] + fak2* S[1][6][3]
                S[3][16][3] = gma1*S[1][16][3] + fak * S[1][16][1]
                S[4][16][3] = gma2*S[1][16][3] + fak * S[1][8][3]
                S[2][16][4] = gma0*S[1][16][4] + fak2* S[1][6][4]
                S[3][16][4] = gma1*S[1][16][4]
                S[4][16][4] = gma2*S[1][16][4] + fak * (S[1][8][4] +S[1][16][1] )
                S[2][17][2] = gma0*S[1][17][2] + fak * (S[1][10][2] +S[1][17][1] )
                S[3][17][2] = gma1*S[1][17][2]
                S[4][17][2] = gma2*S[1][17][2] + fak2* S[1][6][2]
                S[2][17][3] = gma0*S[1][17][3] + fak * S[1][10][3]
                S[3][17][3] = gma1*S[1][17][3] + fak * S[1][17][1]
                S[4][17][3] = gma2*S[1][17][3] + fak2* S[1][6][3]
                S[2][17][4] = gma0*S[1][17][4] + fak * S[1][10][4]
                S[3][17][4] = gma1*S[1][17][4]
                S[4][17][4] = gma2*S[1][17][4] + fak * (2.d0*S[1][6][4] +S[1][17][1] )
                S[2][18][2] = gma0*S[1][18][2] + fak * S[1][18][1]
                S[3][18][2] = gma1*S[1][18][2] + fak2* S[1][7][2]
                S[4][18][2] = gma2*S[1][18][2] + fak * S[1][9][2]
                S[2][18][3] = gma0*S[1][18][3]
                S[3][18][3] = gma1*S[1][18][3] + fak * (2.d0*S[1][7][3] +S[1][18][1] )
                S[4][18][3] = gma2*S[1][18][3] + fak * S[1][9][3]
                S[2][18][4] = gma0*S[1][18][4]
                S[3][18][4] = gma1*S[1][18][4] + fak2* S[1][7][4]
                S[4][18][4] = gma2*S[1][18][4] + fak * (S[1][9][4] +S[1][18][1] )
                S[2][19][2] = gma0*S[1][19][2] + fak * S[1][19][1]
                S[3][19][2] = gma1*S[1][19][2] + fak * S[1][10][2]
                S[4][19][2] = gma2*S[1][19][2] + fak2* S[1][7][2]
                S[2][19][3] = gma0*S[1][19][3]
                S[3][19][3] = gma1*S[1][19][3] + fak * (S[1][10][3] +S[1][19][1] )
                S[4][19][3] = gma2*S[1][19][3] + fak2* S[1][7][3]
                S[2][19][4] = gma0*S[1][19][4]
                S[3][19][4] = gma1*S[1][19][4] + fak * S[1][10][4]
                S[4][19][4] = gma2*S[1][19][4] + fak * (2.d0*S[1][7][4] +S[1][19][1] )
                S[2][20][2] = gma0*S[1][20][2] + fak * (S[1][7][2] +S[1][20][1] )
                S[3][20][2] = gma1*S[1][20][2] + fak * S[1][6][2]
                S[4][20][2] = gma2*S[1][20][2] + fak * S[1][5][2]
                S[2][20][3] = gma0*S[1][20][3] + fak * S[1][7][3]
                S[3][20][3] = gma1*S[1][20][3] + fak * (S[1][6][3] +S[1][20][1] )
                S[4][20][3] = gma2*S[1][20][3] + fak * S[1][5][3]
                S[2][20][4] = gma0*S[1][20][4] + fak * S[1][7][4]
                S[3][20][4] = gma1*S[1][20][4] + fak * S[1][6][4]
                S[4][20][4] = gma2*S[1][20][4] + fak * (S[1][5][4] +S[1][20][1] )
             }

             if ( _lmax_alpha >= 2 && _lmax_gw >= 3 && _lmax_gamma >= 0) {
                S[8][11][1] = gma0*S[2][11][1] + fak * (S[1][11][1] +3.d0*S[2][8][1] )
                S[5][11][1] = gma1*S[2][11][1]
                S[6][11][1] = gma2*S[2][11][1]
                S[8][12][1] = gma0*S[2][12][1] + fak * S[1][12][1]
                S[5][12][1] = gma1*S[2][12][1] + fak3* S[2][9][1]
                S[6][12][1] = gma2*S[2][12][1]
                S[8][13][1] = gma0*S[2][13][1] + fak * S[1][13][1]
                S[5][13][1] = gma1*S[2][13][1]
                S[6][13][1] = gma2*S[2][13][1] + fak3* S[2][10][1]
                S[8][14][1] = gma0*S[2][14][1] + fak * (S[1][14][1] +2.d0*S[2][5][1] )
                S[5][14][1] = gma1*S[2][14][1] + fak * S[2][8][1]
                S[6][14][1] = gma2*S[2][14][1]
                S[8][15][1] = gma0*S[2][15][1] + fak * (S[1][15][1] +S[2][9][1] )
                S[5][15][1] = gma1*S[2][15][1] + fak2* S[2][5][1]
                S[6][15][1] = gma2*S[2][15][1]
                S[8][16][1] = gma0*S[2][16][1] + fak * (S[1][16][1] +2.d0*S[2][6][1] )
                S[5][16][1] = gma1*S[2][16][1]
                S[6][16][1] = gma2*S[2][16][1] + fak * S[2][8][1]
                S[8][17][1] = gma0*S[2][17][1] + fak * (S[1][17][1] +S[2][10][1] )
                S[5][17][1] = gma1*S[2][17][1]
                S[6][17][1] = gma2*S[2][17][1] + fak2* S[2][6][1]
                S[8][18][1] = gma0*S[2][18][1] + fak * S[1][18][1]
                S[5][18][1] = gma1*S[2][18][1] + fak2* S[2][7][1]
                S[6][18][1] = gma2*S[2][18][1] + fak * S[2][9][1]
                S[8][19][1] = gma0*S[2][19][1] + fak * S[1][19][1]
                S[5][19][1] = gma1*S[2][19][1] + fak * S[2][10][1]
                S[6][19][1] = gma2*S[2][19][1] + fak2* S[2][7][1]
                S[8][20][1] = gma0*S[2][20][1] + fak * (S[1][20][1] +S[2][7][1] )
                S[5][20][1] = gma1*S[2][20][1] + fak * S[2][6][1]
                S[6][20][1] = gma2*S[2][20][1] + fak * S[2][5][1]
                S[9][11][1] = gma1*S[3][11][1] + fak * S[1][11][1]
                S[7][11][1] = gma2*S[3][11][1]
                S[9][12][1] = gma1*S[3][12][1] + fak * (S[1][12][1] +3.d0*S[3][9][1] )
                S[7][12][1] = gma2*S[3][12][1]
                S[9][13][1] = gma1*S[3][13][1] + fak * S[1][13][1]
                S[7][13][1] = gma2*S[3][13][1] + fak3* S[3][10][1]
                S[9][14][1] = gma1*S[3][14][1] + fak * (S[1][14][1] +S[3][8][1] )
                S[7][14][1] = gma2*S[3][14][1]
                S[9][15][1] = gma1*S[3][15][1] + fak * (S[1][15][1] +2.d0*S[3][5][1] )
                S[7][15][1] = gma2*S[3][15][1]
                S[9][16][1] = gma1*S[3][16][1] + fak * S[1][16][1]
                S[7][16][1] = gma2*S[3][16][1] + fak * S[3][8][1]
                S[9][17][1] = gma1*S[3][17][1] + fak * S[1][17][1]
                S[7][17][1] = gma2*S[3][17][1] + fak2* S[3][6][1]
                S[9][18][1] = gma1*S[3][18][1] + fak * (S[1][18][1] +2.d0*S[3][7][1] )
                S[7][18][1] = gma2*S[3][18][1] + fak * S[3][9][1]
                S[9][19][1] = gma1*S[3][19][1] + fak * (S[1][19][1] +S[3][10][1] )
                S[7][19][1] = gma2*S[3][19][1] + fak2* S[3][7][1]
                S[9][20][1] = gma1*S[3][20][1] + fak * (S[1][20][1] +S[3][6][1] )
                S[7][20][1] = gma2*S[3][20][1] + fak * S[3][5][1]
                S[10][11][1] = gma2*S[4][11][1] + fak * S[1][11][1]
                S[10][12][1] = gma2*S[4][12][1] + fak * S[1][12][1]
                S[10][13][1] = gma2*S[4][13][1] + fak * (S[1][13][1] +3.d0*S[4][10][1] )
                S[10][14][1] = gma2*S[4][14][1] + fak * S[1][14][1]
                S[10][15][1] = gma2*S[4][15][1] + fak * S[1][15][1]
                S[10][16][1] = gma2*S[4][16][1] + fak * (S[1][16][1] +S[4][8][1] )
                S[10][17][1] = gma2*S[4][17][1] + fak * (S[1][17][1] +2.d0*S[4][6][1] )
                S[10][18][1] = gma2*S[4][18][1] + fak * (S[1][18][1] +S[4][9][1] )
                S[10][19][1] = gma2*S[4][19][1] + fak * (S[1][19][1] +2.d0*S[4][7][1] )
                S[10][20][1] = gma2*S[4][20][1] + fak * (S[1][20][1] +S[4][5][1] )
             }

             if ( _lmax_alpha >= 1 && _lmax_gw >= 3 && _lmax_gamma >= 2) {
                S[2][11][5] = gma0*S[1][11][5] + fak * (3.d0*S[1][8][5] +S[1][11][3] )
                S[3][11][5] = gma1*S[1][11][5] + fak * S[1][11][2]
                S[4][11][5] = gma2*S[1][11][5]
                S[2][11][6] = gma0*S[1][11][6] + fak * (3.d0*S[1][8][6] +S[1][11][4] )
                S[3][11][6] = gma1*S[1][11][6]
                S[4][11][6] = gma2*S[1][11][6] + fak * S[1][11][2]
                S[2][11][7] = gma0*S[1][11][7] + fak3* S[1][8][7]
                S[3][11][7] = gma1*S[1][11][7] + fak * S[1][11][4]
                S[4][11][7] = gma2*S[1][11][7] + fak * S[1][11][3]
                S[2][11][8] = gma0*S[1][11][8] + fak *(3.d0*S[1][8][8]+2.d0*S[1][11][2])
                S[3][11][8] = gma1*S[1][11][8]
                S[4][11][8] = gma2*S[1][11][8]
                S[2][11][9] = gma0*S[1][11][9] + fak3* S[1][8][9]
                S[3][11][9] = gma1*S[1][11][9] + fak2* S[1][11][3]
                S[4][11][9] = gma2*S[1][11][9]
                S[2][11][10] = gma0*S[1][11][10] + fak3* S[1][8][10]
                S[3][11][10] = gma1*S[1][11][10]
                S[4][11][10] = gma2*S[1][11][10] + fak2* S[1][11][4]
                S[2][12][5] = gma0*S[1][12][5] + fak * S[1][12][3]
                S[3][12][5] = gma1*S[1][12][5] + fak * (3.d0*S[1][9][5] +S[1][12][2] )
                S[4][12][5] = gma2*S[1][12][5]
                S[2][12][6] = gma0*S[1][12][6] + fak * S[1][12][4]
                S[3][12][6] = gma1*S[1][12][6] + fak3* S[1][9][6]
                S[4][12][6] = gma2*S[1][12][6] + fak * S[1][12][2]
                S[2][12][7] = gma0*S[1][12][7]
                S[3][12][7] = gma1*S[1][12][7] + fak * (3.d0*S[1][9][7] +S[1][12][4] )
                S[4][12][7] = gma2*S[1][12][7] + fak * S[1][12][3]
                S[2][12][8] = gma0*S[1][12][8] + fak2* S[1][12][2]
                S[3][12][8] = gma1*S[1][12][8] + fak3* S[1][9][8]
                S[4][12][8] = gma2*S[1][12][8]
                S[2][12][9] = gma0*S[1][12][9]
                S[3][12][9] = gma1*S[1][12][9] + fak *(3.d0*S[1][9][9]+2.d0*S[1][12][3])
                S[4][12][9] = gma2*S[1][12][9]
                S[2][12][10] = gma0*S[1][12][10]
                S[3][12][10] = gma1*S[1][12][10] + fak3* S[1][9][10]
                S[4][12][10] = gma2*S[1][12][10] + fak2* S[1][12][4]
                S[2][13][5] = gma0*S[1][13][5] + fak * S[1][13][3]
                S[3][13][5] = gma1*S[1][13][5] + fak * S[1][13][2]
                S[4][13][5] = gma2*S[1][13][5] + fak3* S[1][10][5]
                S[2][13][6] = gma0*S[1][13][6] + fak * S[1][13][4]
                S[3][13][6] = gma1*S[1][13][6]
                S[4][13][6] = gma2*S[1][13][6] + fak * (3.d0*S[1][10][6] +S[1][13][2] )
                S[2][13][7] = gma0*S[1][13][7]
                S[3][13][7] = gma1*S[1][13][7] + fak * S[1][13][4]
                S[4][13][7] = gma2*S[1][13][7] + fak * (3.d0*S[1][10][7] +S[1][13][3] )
                S[2][13][8] = gma0*S[1][13][8] + fak2* S[1][13][2]
                S[3][13][8] = gma1*S[1][13][8]
                S[4][13][8] = gma2*S[1][13][8] + fak3* S[1][10][8]
                S[2][13][9] = gma0*S[1][13][9]
                S[3][13][9] = gma1*S[1][13][9] + fak2* S[1][13][3]
                S[4][13][9] = gma2*S[1][13][9] + fak3* S[1][10][9]
                S[2][13][10] = gma0*S[1][13][10]
                S[3][13][10] = gma1*S[1][13][10]
                S[4][13][10] = gma2*S[1][13][10]  & 
                     +    fak * (3.d0*S[1][10][10] +2.d0*S[1][13][4] )
                S[2][14][5] = gma0*S[1][14][5] + fak * (2.d0*S[1][5][5] +S[1][14][3] )
                S[3][14][5] = gma1*S[1][14][5] + fak * (S[1][8][5] +S[1][14][2] )
                S[4][14][5] = gma2*S[1][14][5]
                S[2][14][6] = gma0*S[1][14][6] + fak * (2.d0*S[1][5][6] +S[1][14][4] )
                S[3][14][6] = gma1*S[1][14][6] + fak * S[1][8][6]
                S[4][14][6] = gma2*S[1][14][6] + fak * S[1][14][2]
                S[2][14][7] = gma0*S[1][14][7] + fak2* S[1][5][7]
                S[3][14][7] = gma1*S[1][14][7] + fak * (S[1][8][7] +S[1][14][4] )
                S[4][14][7] = gma2*S[1][14][7] + fak * S[1][14][3]
                S[2][14][8] = gma0*S[1][14][8] + fak2* (S[1][5][8] +S[1][14][2] )
                S[3][14][8] = gma1*S[1][14][8] + fak * S[1][8][8]
                S[4][14][8] = gma2*S[1][14][8]
                S[2][14][9] = gma0*S[1][14][9] + fak2* S[1][5][9]
                S[3][14][9] = gma1*S[1][14][9] + fak * (S[1][8][9] +2.d0*S[1][14][3] )
                S[4][14][9] = gma2*S[1][14][9]
                S[2][14][10] = gma0*S[1][14][10] + fak2* S[1][5][10]
                S[3][14][10] = gma1*S[1][14][10] + fak * S[1][8][10]
                S[4][14][10] = gma2*S[1][14][10] + fak2* S[1][14][4]
                S[2][15][5] = gma0*S[1][15][5] + fak * (S[1][9][5] +S[1][15][3] )
                S[3][15][5] = gma1*S[1][15][5] + fak * (2.d0*S[1][5][5] +S[1][15][2] )
                S[4][15][5] = gma2*S[1][15][5]
                S[2][15][6] = gma0*S[1][15][6] + fak * (S[1][9][6] +S[1][15][4] )
                S[3][15][6] = gma1*S[1][15][6] + fak2* S[1][5][6]
                S[4][15][6] = gma2*S[1][15][6] + fak * S[1][15][2]
                S[2][15][7] = gma0*S[1][15][7] + fak * S[1][9][7]
                S[3][15][7] = gma1*S[1][15][7] + fak * (2.d0*S[1][5][7] +S[1][15][4] )
                S[4][15][7] = gma2*S[1][15][7] + fak * S[1][15][3]
                S[2][15][8] = gma0*S[1][15][8] + fak * (S[1][9][8] +2.d0*S[1][15][2] )
                S[3][15][8] = gma1*S[1][15][8] + fak2* S[1][5][8]
                S[4][15][8] = gma2*S[1][15][8]
                S[2][15][9] = gma0*S[1][15][9] + fak * S[1][9][9]
                S[3][15][9] = gma1*S[1][15][9] + fak2* (S[1][5][9] +S[1][15][3] )
                S[4][15][9] = gma2*S[1][15][9]
                S[2][15][10] = gma0*S[1][15][10] + fak * S[1][9][10]
                S[3][15][10] = gma1*S[1][15][10] + fak2* S[1][5][10]
                S[4][15][10] = gma2*S[1][15][10] + fak2* S[1][15][4]
                S[2][16][5] = gma0*S[1][16][5] + fak * (2.d0*S[1][6][5] +S[1][16][3] )
                S[3][16][5] = gma1*S[1][16][5] + fak * S[1][16][2]
                S[4][16][5] = gma2*S[1][16][5] + fak * S[1][8][5]
                S[2][16][6] = gma0*S[1][16][6] + fak * (2.d0*S[1][6][6] +S[1][16][4] )
                S[3][16][6] = gma1*S[1][16][6]
                S[4][16][6] = gma2*S[1][16][6] + fak * (S[1][8][6] +S[1][16][2] )
                S[2][16][7] = gma0*S[1][16][7] + fak2* S[1][6][7]
                S[3][16][7] = gma1*S[1][16][7] + fak * S[1][16][4]
                S[4][16][7] = gma2*S[1][16][7] + fak * (S[1][8][7] +S[1][16][3] )
                S[2][16][8] = gma0*S[1][16][8] + fak2* (S[1][6][8] +S[1][16][2] )
                S[3][16][8] = gma1*S[1][16][8]
                S[4][16][8] = gma2*S[1][16][8] + fak * S[1][8][8]
                S[2][16][9] = gma0*S[1][16][9] + fak2* S[1][6][9]
                S[3][16][9] = gma1*S[1][16][9] + fak2* S[1][16][3]
                S[4][16][9] = gma2*S[1][16][9] + fak * S[1][8][9]
                S[2][16][10] = gma0*S[1][16][10] + fak2* S[1][6][10]
                S[3][16][10] = gma1*S[1][16][10]
                S[4][16][10] = gma2*S[1][16][10] + fak * (S[1][8][10] +2.d0*S[1][16][4])
                S[2][17][5] = gma0*S[1][17][5] + fak * (S[1][10][5] +S[1][17][3] )
                S[3][17][5] = gma1*S[1][17][5] + fak * S[1][17][2]
                S[4][17][5] = gma2*S[1][17][5] + fak2* S[1][6][5]
                S[2][17][6] = gma0*S[1][17][6] + fak * (S[1][10][6] +S[1][17][4] )
                S[3][17][6] = gma1*S[1][17][6]
                S[4][17][6] = gma2*S[1][17][6] + fak * (2.d0*S[1][6][6] +S[1][17][2] )
                S[2][17][7] = gma0*S[1][17][7] + fak * S[1][10][7]
                S[3][17][7] = gma1*S[1][17][7] + fak * S[1][17][4]
                S[4][17][7] = gma2*S[1][17][7] + fak * (2.d0*S[1][6][7] +S[1][17][3] )
                S[2][17][8] = gma0*S[1][17][8] + fak * (S[1][10][8] +2.d0*S[1][17][2] )
                S[3][17][8] = gma1*S[1][17][8]
                S[4][17][8] = gma2*S[1][17][8] + fak2* S[1][6][8]
                S[2][17][9] = gma0*S[1][17][9] + fak * S[1][10][9]
                S[3][17][9] = gma1*S[1][17][9] + fak2* S[1][17][3]
                S[4][17][9] = gma2*S[1][17][9] + fak2* S[1][6][9]
                S[2][17][10] = gma0*S[1][17][10] + fak * S[1][10][10]
                S[3][17][10] = gma1*S[1][17][10]
                S[4][17][10] = gma2*S[1][17][10] + fak2* (S[1][6][10] +S[1][17][4] )
                S[2][18][5] = gma0*S[1][18][5] + fak * S[1][18][3]
                S[3][18][5] = gma1*S[1][18][5] + fak * (2.d0*S[1][7][5] +S[1][18][2] )
                S[4][18][5] = gma2*S[1][18][5] + fak * S[1][9][5]
                S[2][18][6] = gma0*S[1][18][6] + fak * S[1][18][4]
                S[3][18][6] = gma1*S[1][18][6] + fak2* S[1][7][6]
                S[4][18][6] = gma2*S[1][18][6] + fak * (S[1][9][6] +S[1][18][2] )
                S[2][18][7] = gma0*S[1][18][7]
                S[3][18][7] = gma1*S[1][18][7] + fak * (2.d0*S[1][7][7] +S[1][18][4] )
                S[4][18][7] = gma2*S[1][18][7] + fak * (S[1][9][7] +S[1][18][3] )
                S[2][18][8] = gma0*S[1][18][8] + fak2* S[1][18][2]
                S[3][18][8] = gma1*S[1][18][8] + fak2* S[1][7][8]
                S[4][18][8] = gma2*S[1][18][8] + fak * S[1][9][8]
                S[2][18][9] = gma0*S[1][18][9]
                S[3][18][9] = gma1*S[1][18][9] + fak2* (S[1][7][9] +S[1][18][3] )
                S[4][18][9] = gma2*S[1][18][9] + fak * S[1][9][9]
                S[2][18][10] = gma0*S[1][18][10]
                S[3][18][10] = gma1*S[1][18][10] + fak2* S[1][7][10]
                S[4][18][10] = gma2*S[1][18][10] + fak * (S[1][9][10] +2.d0*S[1][18][4])
                S[2][19][5] = gma0*S[1][19][5] + fak * S[1][19][3]
                S[3][19][5] = gma1*S[1][19][5] + fak * (S[1][10][5] +S[1][19][2] )
                S[4][19][5] = gma2*S[1][19][5] + fak2* S[1][7][5]
                S[2][19][6] = gma0*S[1][19][6] + fak * S[1][19][4]
                S[3][19][6] = gma1*S[1][19][6] + fak * S[1][10][6]
                S[4][19][6] = gma2*S[1][19][6] + fak * (2.d0*S[1][7][6] +S[1][19][2] )
                S[2][19][7] = gma0*S[1][19][7]
                S[3][19][7] = gma1*S[1][19][7] + fak * (S[1][10][7] +S[1][19][4] )
                S[4][19][7] = gma2*S[1][19][7] + fak * (2.d0*S[1][7][7] +S[1][19][3] )
                S[2][19][8] = gma0*S[1][19][8] + fak2* S[1][19][2]
                S[3][19][8] = gma1*S[1][19][8] + fak * S[1][10][8]
                S[4][19][8] = gma2*S[1][19][8] + fak2* S[1][7][8]
                S[2][19][9] = gma0*S[1][19][9]
                S[3][19][9] = gma1*S[1][19][9] + fak * (S[1][10][9] +2.d0*S[1][19][3] )
                S[4][19][9] = gma2*S[1][19][9] + fak2* S[1][7][9]
                S[2][19][10] = gma0*S[1][19][10]
                S[3][19][10] = gma1*S[1][19][10] + fak * S[1][10][10]
                S[4][19][10] = gma2*S[1][19][10] + fak2* (S[1][7][10] +S[1][19][4] )
                S[2][20][5] = gma0*S[1][20][5] + fak * (S[1][7][5] +S[1][20][3] )
                S[3][20][5] = gma1*S[1][20][5] + fak * (S[1][6][5] +S[1][20][2] )
                S[4][20][5] = gma2*S[1][20][5] + fak * S[1][5][5]
                S[2][20][6] = gma0*S[1][20][6] + fak * (S[1][7][6] +S[1][20][4] )
                S[3][20][6] = gma1*S[1][20][6] + fak * S[1][6][6]
                S[4][20][6] = gma2*S[1][20][6] + fak * (S[1][5][6] +S[1][20][2] )
                S[2][20][7] = gma0*S[1][20][7] + fak * S[1][7][7]
                S[3][20][7] = gma1*S[1][20][7] + fak * (S[1][6][7] +S[1][20][4] )
                S[4][20][7] = gma2*S[1][20][7] + fak * (S[1][5][7] +S[1][20][3] )
                S[2][20][8] = gma0*S[1][20][8] + fak * (S[1][7][8] +2.d0*S[1][20][2] )
                S[3][20][8] = gma1*S[1][20][8] + fak * S[1][6][8]
                S[4][20][8] = gma2*S[1][20][8] + fak * S[1][5][8]
                S[2][20][9] = gma0*S[1][20][9] + fak * S[1][7][9]
                S[3][20][9] = gma1*S[1][20][9] + fak * (S[1][6][9] +2.d0*S[1][20][3] )
                S[4][20][9] = gma2*S[1][20][9] + fak * S[1][5][9]
                S[2][20][10] = gma0*S[1][20][10] + fak * S[1][7][10]
                S[3][20][10] = gma1*S[1][20][10] + fak * S[1][6][10]
                S[4][20][10] = gma2*S[1][20][10] + fak * (S[1][5][10] +2.d0*S[1][20][4])
             }

             if ( _lmax_alpha >= 2 && _lmax_gw >= 3 && _lmax_gamma >= 1) {
                S[8][11][2] = gma0*S[2][11][2]  & 
                     +    fak * (S[1][11][2] +3.d0*S[2][8][2] +S[2][11][1] )
                S[5][11][2] = gma1*S[2][11][2]
                S[6][11][2] = gma2*S[2][11][2]
                S[8][11][3] = gma0*S[2][11][3] + fak * (S[1][11][3] +3.d0*S[2][8][3] )
                S[5][11][3] = gma1*S[2][11][3] + fak * S[2][11][1]
                S[6][11][3] = gma2*S[2][11][3]
                S[8][11][4] = gma0*S[2][11][4] + fak * (S[1][11][4] +3.d0*S[2][8][4] )
                S[5][11][4] = gma1*S[2][11][4]
                S[6][11][4] = gma2*S[2][11][4] + fak * S[2][11][1]
                S[8][12][2] = gma0*S[2][12][2] + fak * (S[1][12][2] +S[2][12][1] )
                S[5][12][2] = gma1*S[2][12][2] + fak3* S[2][9][2]
                S[6][12][2] = gma2*S[2][12][2]
                S[8][12][3] = gma0*S[2][12][3] + fak * S[1][12][3]
                S[5][12][3] = gma1*S[2][12][3] + fak * (3.d0*S[2][9][3] +S[2][12][1] )
                S[6][12][3] = gma2*S[2][12][3]
                S[8][12][4] = gma0*S[2][12][4] + fak * S[1][12][4]
                S[5][12][4] = gma1*S[2][12][4] + fak3* S[2][9][4]
                S[6][12][4] = gma2*S[2][12][4] + fak * S[2][12][1]
                S[8][13][2] = gma0*S[2][13][2] + fak * (S[1][13][2] +S[2][13][1] )
                S[5][13][2] = gma1*S[2][13][2]
                S[6][13][2] = gma2*S[2][13][2] + fak3* S[2][10][2]
                S[8][13][3] = gma0*S[2][13][3] + fak * S[1][13][3]
                S[5][13][3] = gma1*S[2][13][3] + fak * S[2][13][1]
                S[6][13][3] = gma2*S[2][13][3] + fak3* S[2][10][3]
                S[8][13][4] = gma0*S[2][13][4] + fak * S[1][13][4]
                S[5][13][4] = gma1*S[2][13][4]
                S[6][13][4] = gma2*S[2][13][4] + fak * (3.d0*S[2][10][4] +S[2][13][1] )
                S[8][14][2] = gma0*S[2][14][2]  & 
                     +    fak * (S[1][14][2] +2.d0*S[2][5][2] +S[2][14][1] )
                S[5][14][2] = gma1*S[2][14][2] + fak * S[2][8][2]
                S[6][14][2] = gma2*S[2][14][2]
                S[8][14][3] = gma0*S[2][14][3] + fak * (S[1][14][3] +2.d0*S[2][5][3] )
                S[5][14][3] = gma1*S[2][14][3] + fak * (S[2][8][3] +S[2][14][1] )
                S[6][14][3] = gma2*S[2][14][3]
                S[8][14][4] = gma0*S[2][14][4] + fak * (S[1][14][4] +2.d0*S[2][5][4] )
                S[5][14][4] = gma1*S[2][14][4] + fak * S[2][8][4]
                S[6][14][4] = gma2*S[2][14][4] + fak * S[2][14][1]
                S[8][15][2] = gma0*S[2][15][2] + fak *(S[1][15][2]+S[2][9][2]+S[2][15][1])
                S[5][15][2] = gma1*S[2][15][2] + fak2* S[2][5][2]
                S[6][15][2] = gma2*S[2][15][2]
                S[8][15][3] = gma0*S[2][15][3] + fak * (S[1][15][3] +S[2][9][3] )
                S[5][15][3] = gma1*S[2][15][3] + fak * (2.d0*S[2][5][3] +S[2][15][1] )
                S[6][15][3] = gma2*S[2][15][3]
                S[8][15][4] = gma0*S[2][15][4] + fak * (S[1][15][4] +S[2][9][4] )
                S[5][15][4] = gma1*S[2][15][4] + fak2* S[2][5][4]
                S[6][15][4] = gma2*S[2][15][4] + fak * S[2][15][1]
                S[8][16][2] = gma0*S[2][16][2]  & 
                     +    fak * (S[1][16][2] +2.d0*S[2][6][2] +S[2][16][1] )
                S[5][16][2] = gma1*S[2][16][2]
                S[6][16][2] = gma2*S[2][16][2] + fak * S[2][8][2]
                S[8][16][3] = gma0*S[2][16][3] + fak * (S[1][16][3] +2.d0*S[2][6][3] )
                S[5][16][3] = gma1*S[2][16][3] + fak * S[2][16][1]
                S[6][16][3] = gma2*S[2][16][3] + fak * S[2][8][3]
                S[8][16][4] = gma0*S[2][16][4] + fak * (S[1][16][4] +2.d0*S[2][6][4] )
                S[5][16][4] = gma1*S[2][16][4]
                S[6][16][4] = gma2*S[2][16][4] + fak * (S[2][8][4] +S[2][16][1] )
                S[8][17][2] = gma0*S[2][17][2] + fak*(S[1][17][2]+S[2][10][2]+S[2][17][1])
                S[5][17][2] = gma1*S[2][17][2]
                S[6][17][2] = gma2*S[2][17][2] + fak2* S[2][6][2]
                S[8][17][3] = gma0*S[2][17][3] + fak * (S[1][17][3] +S[2][10][3] )
                S[5][17][3] = gma1*S[2][17][3] + fak * S[2][17][1]
                S[6][17][3] = gma2*S[2][17][3] + fak2* S[2][6][3]
                S[8][17][4] = gma0*S[2][17][4] + fak * (S[1][17][4] +S[2][10][4] )
                S[5][17][4] = gma1*S[2][17][4]
                S[6][17][4] = gma2*S[2][17][4] + fak * (2.d0*S[2][6][4] +S[2][17][1] )
                S[8][18][2] = gma0*S[2][18][2] + fak * (S[1][18][2] +S[2][18][1] )
                S[5][18][2] = gma1*S[2][18][2] + fak2* S[2][7][2]
                S[6][18][2] = gma2*S[2][18][2] + fak * S[2][9][2]
                S[8][18][3] = gma0*S[2][18][3] + fak * S[1][18][3]
                S[5][18][3] = gma1*S[2][18][3] + fak * (2.d0*S[2][7][3] +S[2][18][1] )
                S[6][18][3] = gma2*S[2][18][3] + fak * S[2][9][3]
                S[8][18][4] = gma0*S[2][18][4] + fak * S[1][18][4]
                S[5][18][4] = gma1*S[2][18][4] + fak2* S[2][7][4]
                S[6][18][4] = gma2*S[2][18][4] + fak * (S[2][9][4] +S[2][18][1] )
                S[8][19][2] = gma0*S[2][19][2] + fak * (S[1][19][2] +S[2][19][1] )
                S[5][19][2] = gma1*S[2][19][2] + fak * S[2][10][2]
                S[6][19][2] = gma2*S[2][19][2] + fak2* S[2][7][2]
                S[8][19][3] = gma0*S[2][19][3] + fak * S[1][19][3]
                S[5][19][3] = gma1*S[2][19][3] + fak * (S[2][10][3] +S[2][19][1] )
                S[6][19][3] = gma2*S[2][19][3] + fak2* S[2][7][3]
                S[8][19][4] = gma0*S[2][19][4] + fak * S[1][19][4]
                S[5][19][4] = gma1*S[2][19][4] + fak * S[2][10][4]
                S[6][19][4] = gma2*S[2][19][4] + fak * (2.d0*S[2][7][4] +S[2][19][1] )
                S[8][20][2] = gma0*S[2][20][2] + fak *(S[1][20][2]+S[2][7][2]+S[2][20][1])
                S[5][20][2] = gma1*S[2][20][2] + fak * S[2][6][2]
                S[6][20][2] = gma2*S[2][20][2] + fak * S[2][5][2]
                S[8][20][3] = gma0*S[2][20][3] + fak * (S[1][20][3] +S[2][7][3] )
                S[5][20][3] = gma1*S[2][20][3] + fak * (S[2][6][3] +S[2][20][1] )
                S[6][20][3] = gma2*S[2][20][3] + fak * S[2][5][3]
                S[8][20][4] = gma0*S[2][20][4] + fak * (S[1][20][4] +S[2][7][4] )
                S[5][20][4] = gma1*S[2][20][4] + fak * S[2][6][4]
                S[6][20][4] = gma2*S[2][20][4] + fak * (S[2][5][4] +S[2][20][1] )
                S[9][11][2] = gma1*S[3][11][2] + fak * S[1][11][2]
                S[7][11][2] = gma2*S[3][11][2]
                S[9][11][3] = gma1*S[3][11][3] + fak * (S[1][11][3] +S[3][11][1] )
                S[7][11][3] = gma2*S[3][11][3]
                S[9][11][4] = gma1*S[3][11][4] + fak * S[1][11][4]
                S[7][11][4] = gma2*S[3][11][4] + fak * S[3][11][1]
                S[9][12][2] = gma1*S[3][12][2] + fak * (S[1][12][2] +3.d0*S[3][9][2] )
                S[7][12][2] = gma2*S[3][12][2]
                S[9][12][3] = gma1*S[3][12][3]  & 
                     +    fak * (S[1][12][3] +3.d0*S[3][9][3] +S[3][12][1] )
                S[7][12][3] = gma2*S[3][12][3]
                S[9][12][4] = gma1*S[3][12][4] + fak * (S[1][12][4] +3.d0*S[3][9][4] )
                S[7][12][4] = gma2*S[3][12][4] + fak * S[3][12][1]
                S[9][13][2] = gma1*S[3][13][2] + fak * S[1][13][2]
                S[7][13][2] = gma2*S[3][13][2] + fak3* S[3][10][2]
                S[9][13][3] = gma1*S[3][13][3] + fak * (S[1][13][3] +S[3][13][1] )
                S[7][13][3] = gma2*S[3][13][3] + fak3* S[3][10][3]
                S[9][13][4] = gma1*S[3][13][4] + fak * S[1][13][4]
                S[7][13][4] = gma2*S[3][13][4] + fak * (3.d0*S[3][10][4] +S[3][13][1] )
                S[9][14][2] = gma1*S[3][14][2] + fak * (S[1][14][2] +S[3][8][2] )
                S[7][14][2] = gma2*S[3][14][2]
                S[9][14][3] = gma1*S[3][14][3] + fak *(S[1][14][3]+S[3][8][3]+S[3][14][1])
                S[7][14][3] = gma2*S[3][14][3]
                S[9][14][4] = gma1*S[3][14][4] + fak * (S[1][14][4] +S[3][8][4] )
                S[7][14][4] = gma2*S[3][14][4] + fak * S[3][14][1]
                S[9][15][2] = gma1*S[3][15][2] + fak * (S[1][15][2] +2.d0*S[3][5][2] )
                S[7][15][2] = gma2*S[3][15][2]
                S[9][15][3] = gma1*S[3][15][3]  & 
                     +    fak * (S[1][15][3] +2.d0*S[3][5][3] +S[3][15][1] )
                S[7][15][3] = gma2*S[3][15][3]
                S[9][15][4] = gma1*S[3][15][4] + fak * (S[1][15][4] +2.d0*S[3][5][4] )
                S[7][15][4] = gma2*S[3][15][4] + fak * S[3][15][1]
                S[9][16][2] = gma1*S[3][16][2] + fak * S[1][16][2]
                S[7][16][2] = gma2*S[3][16][2] + fak * S[3][8][2]
                S[9][16][3] = gma1*S[3][16][3] + fak * (S[1][16][3] +S[3][16][1] )
                S[7][16][3] = gma2*S[3][16][3] + fak * S[3][8][3]
                S[9][16][4] = gma1*S[3][16][4] + fak * S[1][16][4]
                S[7][16][4] = gma2*S[3][16][4] + fak * (S[3][8][4] +S[3][16][1] )
                S[9][17][2] = gma1*S[3][17][2] + fak * S[1][17][2]
                S[7][17][2] = gma2*S[3][17][2] + fak2* S[3][6][2]
                S[9][17][3] = gma1*S[3][17][3] + fak * (S[1][17][3] +S[3][17][1] )
                S[7][17][3] = gma2*S[3][17][3] + fak2* S[3][6][3]
                S[9][17][4] = gma1*S[3][17][4] + fak * S[1][17][4]
                S[7][17][4] = gma2*S[3][17][4] + fak * (2.d0*S[3][6][4] +S[3][17][1] )
                S[9][18][2] = gma1*S[3][18][2] + fak * (S[1][18][2] +2.d0*S[3][7][2] )
                S[7][18][2] = gma2*S[3][18][2] + fak * S[3][9][2]
                S[9][18][3] = gma1*S[3][18][3]  & 
                     +    fak * (S[1][18][3] +2.d0*S[3][7][3] +S[3][18][1] )
                S[7][18][3] = gma2*S[3][18][3] + fak * S[3][9][3]
                S[9][18][4] = gma1*S[3][18][4] + fak * (S[1][18][4] +2.d0*S[3][7][4] )
                S[7][18][4] = gma2*S[3][18][4] + fak * (S[3][9][4] +S[3][18][1] )
                S[9][19][2] = gma1*S[3][19][2] + fak * (S[1][19][2] +S[3][10][2] )
                S[7][19][2] = gma2*S[3][19][2] + fak2* S[3][7][2]
                S[9][19][3] = gma1*S[3][19][3] + fak*(S[1][19][3]+S[3][10][3]+S[3][19][1])
                S[7][19][3] = gma2*S[3][19][3] + fak2* S[3][7][3]
                S[9][19][4] = gma1*S[3][19][4] + fak * (S[1][19][4] +S[3][10][4] )
                S[7][19][4] = gma2*S[3][19][4] + fak * (2.d0*S[3][7][4] +S[3][19][1] )
                S[9][20][2] = gma1*S[3][20][2] + fak * (S[1][20][2] +S[3][6][2] )
                S[7][20][2] = gma2*S[3][20][2] + fak * S[3][5][2]
                S[9][20][3] = gma1*S[3][20][3] + fak *(S[1][20][3]+S[3][6][3]+S[3][20][1])
                S[7][20][3] = gma2*S[3][20][3] + fak * S[3][5][3]
                S[9][20][4] = gma1*S[3][20][4] + fak * (S[1][20][4] +S[3][6][4] )
                S[7][20][4] = gma2*S[3][20][4] + fak * (S[3][5][4] +S[3][20][1] )
                S[10][11][2] = gma2*S[4][11][2] + fak * S[1][11][2]
                S[10][11][3] = gma2*S[4][11][3] + fak * S[1][11][3]
                S[10][11][4] = gma2*S[4][11][4] + fak * (S[1][11][4] +S[4][11][1] )
                S[10][12][2] = gma2*S[4][12][2] + fak * S[1][12][2]
                S[10][12][3] = gma2*S[4][12][3] + fak * S[1][12][3]
                S[10][12][4] = gma2*S[4][12][4] + fak * (S[1][12][4] +S[4][12][1] )
                S[10][13][2] = gma2*S[4][13][2] + fak * (S[1][13][2] +3.d0*S[4][10][2] )
                S[10][13][3] = gma2*S[4][13][3] + fak * (S[1][13][3] +3.d0*S[4][10][3] )
                S[10][13][4] = gma2*S[4][13][4]  & 
                     +    fak * (S[1][13][4] +3.d0*S[4][10][4] +S[4][13][1] )
                S[10][14][2] = gma2*S[4][14][2] + fak * S[1][14][2]
                S[10][14][3] = gma2*S[4][14][3] + fak * S[1][14][3]
                S[10][14][4] = gma2*S[4][14][4] + fak * (S[1][14][4] +S[4][14][1] )
                S[10][15][2] = gma2*S[4][15][2] + fak * S[1][15][2]
                S[10][15][3] = gma2*S[4][15][3] + fak * S[1][15][3]
                S[10][15][4] = gma2*S[4][15][4] + fak * (S[1][15][4] +S[4][15][1] )
                S[10][16][2] = gma2*S[4][16][2] + fak * (S[1][16][2] +S[4][8][2] )
                S[10][16][3] = gma2*S[4][16][3] + fak * (S[1][16][3] +S[4][8][3] )
                S[10][16][4] = gma2*S[4][16][4] + fak*(S[1][16][4]+S[4][8][4]+S[4][16][1])
                S[10][17][2] = gma2*S[4][17][2] + fak * (S[1][17][2] +2.d0*S[4][6][2] )
                S[10][17][3] = gma2*S[4][17][3] + fak * (S[1][17][3] +2.d0*S[4][6][3] )
                S[10][17][4] = gma2*S[4][17][4]  & 
                     +    fak * (S[1][17][4] +2.d0*S[4][6][4] +S[4][17][1] )
                S[10][18][2] = gma2*S[4][18][2] + fak * (S[1][18][2] +S[4][9][2] )
                S[10][18][3] = gma2*S[4][18][3] + fak * (S[1][18][3] +S[4][9][3] )
                S[10][18][4] = gma2*S[4][18][4] + fak*(S[1][18][4]+S[4][9][4]+S[4][18][1])
                S[10][19][2] = gma2*S[4][19][2] + fak * (S[1][19][2] +2.d0*S[4][7][2] )
                S[10][19][3] = gma2*S[4][19][3] + fak * (S[1][19][3] +2.d0*S[4][7][3] )
                S[10][19][4] = gma2*S[4][19][4]  & 
                     +    fak * (S[1][19][4] +2.d0*S[4][7][4] +S[4][19][1] )
                S[10][20][2] = gma2*S[4][20][2] + fak * (S[1][20][2] +S[4][5][2] )
                S[10][20][3] = gma2*S[4][20][3] + fak * (S[1][20][3] +S[4][5][3] )
                S[10][20][4] = gma2*S[4][20][4] + fak*(S[1][20][4]+S[4][5][4]+S[4][20][1])
             }

             if ( _lmax_alpha >= 2 && _lmax_gw >= 3 && _lmax_gamma >= 2) {
                S[8][11][5] = gma0*S[2][11][5]  & 
                     +    fak * (S[1][11][5] +3.d0*S[2][8][5] +S[2][11][3] )
                S[5][11][5] = gma1*S[2][11][5] + fak * S[2][11][2]
                S[6][11][5] = gma2*S[2][11][5]
                S[8][11][6] = gma0*S[2][11][6]  & 
                     +    fak * (S[1][11][6] +3.d0*S[2][8][6] +S[2][11][4] )
                S[5][11][6] = gma1*S[2][11][6]
                S[6][11][6] = gma2*S[2][11][6] + fak * S[2][11][2]
                S[8][11][7] = gma0*S[2][11][7] + fak * (S[1][11][7] +3.d0*S[2][8][7] )
                S[5][11][7] = gma1*S[2][11][7] + fak * S[2][11][4]
                S[6][11][7] = gma2*S[2][11][7] + fak * S[2][11][3]
                S[8][11][8] = gma0*S[2][11][8]  & 
                     +    fak * (S[1][11][8] +3.d0*S[2][8][8] +2.d0*S[2][11][2] )
                S[5][11][8] = gma1*S[2][11][8]
                S[6][11][8] = gma2*S[2][11][8]
                S[8][11][9] = gma0*S[2][11][9] + fak * (S[1][11][9] +3.d0*S[2][8][9] )
                S[5][11][9] = gma1*S[2][11][9] + fak2* S[2][11][3]
                S[6][11][9] = gma2*S[2][11][9]
                S[8][11][10] = gma0*S[2][11][10] + fak * (S[1][11][10]+3.d0*S[2][8][10])
                S[5][11][10] = gma1*S[2][11][10]
                S[6][11][10] = gma2*S[2][11][10] + fak2* S[2][11][4]
                S[8][12][5] = gma0*S[2][12][5] + fak * (S[1][12][5] +S[2][12][3] )
                S[5][12][5] = gma1*S[2][12][5] + fak * (3.d0*S[2][9][5] +S[2][12][2] )
                S[6][12][5] = gma2*S[2][12][5]
                S[8][12][6] = gma0*S[2][12][6] + fak * (S[1][12][6] +S[2][12][4] )
                S[5][12][6] = gma1*S[2][12][6] + fak3* S[2][9][6]
                S[6][12][6] = gma2*S[2][12][6] + fak * S[2][12][2]
                S[8][12][7] = gma0*S[2][12][7] + fak * S[1][12][7]
                S[5][12][7] = gma1*S[2][12][7] + fak * (3.d0*S[2][9][7] +S[2][12][4] )
                S[6][12][7] = gma2*S[2][12][7] + fak * S[2][12][3]
                S[8][12][8] = gma0*S[2][12][8] + fak * (S[1][12][8] +2.d0*S[2][12][2] )
                S[5][12][8] = gma1*S[2][12][8] + fak3* S[2][9][8]
                S[6][12][8] = gma2*S[2][12][8]
                S[8][12][9] = gma0*S[2][12][9] + fak * S[1][12][9]
                S[5][12][9] = gma1*S[2][12][9] + fak *(3.d0*S[2][9][9]+2.d0*S[2][12][3])
                S[6][12][9] = gma2*S[2][12][9]
                S[8][12][10] = gma0*S[2][12][10] + fak * S[1][12][10]
                S[5][12][10] = gma1*S[2][12][10] + fak3* S[2][9][10]
                S[6][12][10] = gma2*S[2][12][10] + fak2* S[2][12][4]
                S[8][13][5] = gma0*S[2][13][5] + fak * (S[1][13][5] +S[2][13][3] )
                S[5][13][5] = gma1*S[2][13][5] + fak * S[2][13][2]
                S[6][13][5] = gma2*S[2][13][5] + fak3* S[2][10][5]
                S[8][13][6] = gma0*S[2][13][6] + fak * (S[1][13][6] +S[2][13][4] )
                S[5][13][6] = gma1*S[2][13][6]
                S[6][13][6] = gma2*S[2][13][6] + fak * (3.d0*S[2][10][6] +S[2][13][2] )
                S[8][13][7] = gma0*S[2][13][7] + fak * S[1][13][7]
                S[5][13][7] = gma1*S[2][13][7] + fak * S[2][13][4]
                S[6][13][7] = gma2*S[2][13][7] + fak * (3.d0*S[2][10][7] +S[2][13][3] )
                S[8][13][8] = gma0*S[2][13][8] + fak * (S[1][13][8] +2.d0*S[2][13][2] )
                S[5][13][8] = gma1*S[2][13][8]
                S[6][13][8] = gma2*S[2][13][8] + fak3* S[2][10][8]
                S[8][13][9] = gma0*S[2][13][9] + fak * S[1][13][9]
                S[5][13][9] = gma1*S[2][13][9] + fak2* S[2][13][3]
                S[6][13][9] = gma2*S[2][13][9] + fak3* S[2][10][9]
                S[8][13][10] = gma0*S[2][13][10] + fak * S[1][13][10]
                S[5][13][10] = gma1*S[2][13][10]
                S[6][13][10] = gma2*S[2][13][10]  & 
                     +    fak * (3.d0*S[2][10][10] +2.d0*S[2][13][4] )
                S[8][14][5] = gma0*S[2][14][5]  & 
                     +    fak * (S[1][14][5] +2.d0*S[2][5][5] +S[2][14][3] )
                S[5][14][5] = gma1*S[2][14][5] + fak * (S[2][8][5] +S[2][14][2] )
                S[6][14][5] = gma2*S[2][14][5]
                S[8][14][6] = gma0*S[2][14][6]  & 
                     +    fak * (S[1][14][6] +2.d0*S[2][5][6] +S[2][14][4] )
                S[5][14][6] = gma1*S[2][14][6] + fak * S[2][8][6]
                S[6][14][6] = gma2*S[2][14][6] + fak * S[2][14][2]
                S[8][14][7] = gma0*S[2][14][7] + fak * (S[1][14][7] +2.d0*S[2][5][7] )
                S[5][14][7] = gma1*S[2][14][7] + fak * (S[2][8][7] +S[2][14][4] )
                S[6][14][7] = gma2*S[2][14][7] + fak * S[2][14][3]
                S[8][14][8] = gma0*S[2][14][8]  & 
                     +    fak * (S[1][14][8] +2.d0*S[2][5][8] +2.d0*S[2][14][2] )
                S[5][14][8] = gma1*S[2][14][8] + fak * S[2][8][8]
                S[6][14][8] = gma2*S[2][14][8]
                S[8][14][9] = gma0*S[2][14][9] + fak * (S[1][14][9] +2.d0*S[2][5][9] )
                S[5][14][9] = gma1*S[2][14][9] + fak * (S[2][8][9] +2.d0*S[2][14][3] )
                S[6][14][9] = gma2*S[2][14][9]
                S[8][14][10] = gma0*S[2][14][10] + fak * (S[1][14][10]+2.d0*S[2][5][10])
                S[5][14][10] = gma1*S[2][14][10] + fak * S[2][8][10]
                S[6][14][10] = gma2*S[2][14][10] + fak2* S[2][14][4]
                S[8][15][5] = gma0*S[2][15][5] + fak *(S[1][15][5]+S[2][9][5]+S[2][15][3])
                S[5][15][5] = gma1*S[2][15][5] + fak * (2.d0*S[2][5][5] +S[2][15][2] )
                S[6][15][5] = gma2*S[2][15][5]
                S[8][15][6] = gma0*S[2][15][6] + fak *(S[1][15][6]+S[2][9][6]+S[2][15][4])
                S[5][15][6] = gma1*S[2][15][6] + fak2* S[2][5][6]
                S[6][15][6] = gma2*S[2][15][6] + fak * S[2][15][2]
                S[8][15][7] = gma0*S[2][15][7] + fak * (S[1][15][7] +S[2][9][7] )
                S[5][15][7] = gma1*S[2][15][7] + fak * (2.d0*S[2][5][7] +S[2][15][4] )
                S[6][15][7] = gma2*S[2][15][7] + fak * S[2][15][3]
                S[8][15][8] = gma0*S[2][15][8]  & 
                     +    fak * (S[1][15][8] +S[2][9][8] +2.d0*S[2][15][2] )
                S[5][15][8] = gma1*S[2][15][8] + fak2* S[2][5][8]
                S[6][15][8] = gma2*S[2][15][8]
                S[8][15][9] = gma0*S[2][15][9] + fak * (S[1][15][9] +S[2][9][9] )
                S[5][15][9] = gma1*S[2][15][9] + fak2* (S[2][5][9] +S[2][15][3] )
                S[6][15][9] = gma2*S[2][15][9]
                S[8][15][10] = gma0*S[2][15][10] + fak * (S[1][15][10] +S[2][9][10] )
                S[5][15][10] = gma1*S[2][15][10] + fak2* S[2][5][10]
                S[6][15][10] = gma2*S[2][15][10] + fak2* S[2][15][4]
                S[8][16][5] = gma0*S[2][16][5]  & 
                     +    fak * (S[1][16][5] +2.d0*S[2][6][5] +S[2][16][3] )
                S[5][16][5] = gma1*S[2][16][5] + fak * S[2][16][2]
                S[6][16][5] = gma2*S[2][16][5] + fak * S[2][8][5]
                S[8][16][6] = gma0*S[2][16][6]  & 
                     +    fak * (S[1][16][6] +2.d0*S[2][6][6] +S[2][16][4] )
                S[5][16][6] = gma1*S[2][16][6]
                S[6][16][6] = gma2*S[2][16][6] + fak * (S[2][8][6] +S[2][16][2] )
                S[8][16][7] = gma0*S[2][16][7] + fak * (S[1][16][7] +2.d0*S[2][6][7] )
                S[5][16][7] = gma1*S[2][16][7] + fak * S[2][16][4]
                S[6][16][7] = gma2*S[2][16][7] + fak * (S[2][8][7] +S[2][16][3] )
                S[8][16][8] = gma0*S[2][16][8]  & 
                     +    fak * (S[1][16][8] +2.d0*S[2][6][8] +2.d0*S[2][16][2] )
                S[5][16][8] = gma1*S[2][16][8]
                S[6][16][8] = gma2*S[2][16][8] + fak * S[2][8][8]
                S[8][16][9] = gma0*S[2][16][9] + fak * (S[1][16][9] +2.d0*S[2][6][9] )
                S[5][16][9] = gma1*S[2][16][9] + fak2* S[2][16][3]
                S[6][16][9] = gma2*S[2][16][9] + fak * S[2][8][9]
                S[8][16][10] = gma0*S[2][16][10] + fak * (S[1][16][10]+2.d0*S[2][6][10])
                S[5][16][10] = gma1*S[2][16][10]
                S[6][16][10] = gma2*S[2][16][10] + fak * (S[2][8][10] +2.d0*S[2][16][4])
                S[8][17][5] = gma0*S[2][17][5] + fak*(S[1][17][5]+S[2][10][5]+S[2][17][3])
                S[5][17][5] = gma1*S[2][17][5] + fak * S[2][17][2]
                S[6][17][5] = gma2*S[2][17][5] + fak2* S[2][6][5]
                S[8][17][6] = gma0*S[2][17][6] + fak*(S[1][17][6]+S[2][10][6]+S[2][17][4])
                S[5][17][6] = gma1*S[2][17][6]
                S[6][17][6] = gma2*S[2][17][6] + fak * (2.d0*S[2][6][6] +S[2][17][2] )
                S[8][17][7] = gma0*S[2][17][7] + fak * (S[1][17][7] +S[2][10][7] )
                S[5][17][7] = gma1*S[2][17][7] + fak * S[2][17][4]
                S[6][17][7] = gma2*S[2][17][7] + fak * (2.d0*S[2][6][7] +S[2][17][3] )
                S[8][17][8] = gma0*S[2][17][8]  & 
                     +    fak * (S[1][17][8] +S[2][10][8] +2.d0*S[2][17][2] )
                S[5][17][8] = gma1*S[2][17][8]
                S[6][17][8] = gma2*S[2][17][8] + fak2* S[2][6][8]
                S[8][17][9] = gma0*S[2][17][9] + fak * (S[1][17][9] +S[2][10][9] )
                S[5][17][9] = gma1*S[2][17][9] + fak2* S[2][17][3]
                S[6][17][9] = gma2*S[2][17][9] + fak2* S[2][6][9]
                S[8][17][10] = gma0*S[2][17][10] + fak * (S[1][17][10] +S[2][10][10] )
                S[5][17][10] = gma1*S[2][17][10]
                S[6][17][10] = gma2*S[2][17][10] + fak2* (S[2][6][10] +S[2][17][4] )
                S[8][18][5] = gma0*S[2][18][5] + fak * (S[1][18][5] +S[2][18][3] )
                S[5][18][5] = gma1*S[2][18][5] + fak * (2.d0*S[2][7][5] +S[2][18][2] )
                S[6][18][5] = gma2*S[2][18][5] + fak * S[2][9][5]
                S[8][18][6] = gma0*S[2][18][6] + fak * (S[1][18][6] +S[2][18][4] )
                S[5][18][6] = gma1*S[2][18][6] + fak2* S[2][7][6]
                S[6][18][6] = gma2*S[2][18][6] + fak * (S[2][9][6] +S[2][18][2] )
                S[8][18][7] = gma0*S[2][18][7] + fak * S[1][18][7]
                S[5][18][7] = gma1*S[2][18][7] + fak * (2.d0*S[2][7][7] +S[2][18][4] )
                S[6][18][7] = gma2*S[2][18][7] + fak * (S[2][9][7] +S[2][18][3] )
                S[8][18][8] = gma0*S[2][18][8] + fak * (S[1][18][8] +2.d0*S[2][18][2] )
                S[5][18][8] = gma1*S[2][18][8] + fak2* S[2][7][8]
                S[6][18][8] = gma2*S[2][18][8] + fak * S[2][9][8]
                S[8][18][9] = gma0*S[2][18][9] + fak * S[1][18][9]
                S[5][18][9] = gma1*S[2][18][9] + fak2* (S[2][7][9] +S[2][18][3] )
                S[6][18][9] = gma2*S[2][18][9] + fak * S[2][9][9]
                S[8][18][10] = gma0*S[2][18][10] + fak * S[1][18][10]
                S[5][18][10] = gma1*S[2][18][10] + fak2* S[2][7][10]
                S[6][18][10] = gma2*S[2][18][10] + fak * (S[2][9][10] +2.d0*S[2][18][4])
                S[8][19][5] = gma0*S[2][19][5] + fak * (S[1][19][5] +S[2][19][3] )
                S[5][19][5] = gma1*S[2][19][5] + fak * (S[2][10][5] +S[2][19][2] )
                S[6][19][5] = gma2*S[2][19][5] + fak2* S[2][7][5]
                S[8][19][6] = gma0*S[2][19][6] + fak * (S[1][19][6] +S[2][19][4] )
                S[5][19][6] = gma1*S[2][19][6] + fak * S[2][10][6]
                S[6][19][6] = gma2*S[2][19][6] + fak * (2.d0*S[2][7][6] +S[2][19][2] )
                S[8][19][7] = gma0*S[2][19][7] + fak * S[1][19][7]
                S[5][19][7] = gma1*S[2][19][7] + fak * (S[2][10][7] +S[2][19][4] )
                S[6][19][7] = gma2*S[2][19][7] + fak * (2.d0*S[2][7][7] +S[2][19][3] )
                S[8][19][8] = gma0*S[2][19][8] + fak * (S[1][19][8] +2.d0*S[2][19][2] )
                S[5][19][8] = gma1*S[2][19][8] + fak * S[2][10][8]
                S[6][19][8] = gma2*S[2][19][8] + fak2* S[2][7][8]
                S[8][19][9] = gma0*S[2][19][9] + fak * S[1][19][9]
                S[5][19][9] = gma1*S[2][19][9] + fak * (S[2][10][9] +2.d0*S[2][19][3] )
                S[6][19][9] = gma2*S[2][19][9] + fak2* S[2][7][9]
                S[8][19][10] = gma0*S[2][19][10] + fak * S[1][19][10]
                S[5][19][10] = gma1*S[2][19][10] + fak * S[2][10][10]
                S[6][19][10] = gma2*S[2][19][10] + fak2* (S[2][7][10] +S[2][19][4] )
                S[8][20][5] = gma0*S[2][20][5] + fak *(S[1][20][5]+S[2][7][5]+S[2][20][3])
                S[5][20][5] = gma1*S[2][20][5] + fak * (S[2][6][5] +S[2][20][2] )
                S[6][20][5] = gma2*S[2][20][5] + fak * S[2][5][5]
                S[8][20][6] = gma0*S[2][20][6] + fak *(S[1][20][6]+S[2][7][6]+S[2][20][4])
                S[5][20][6] = gma1*S[2][20][6] + fak * S[2][6][6]
                S[6][20][6] = gma2*S[2][20][6] + fak * (S[2][5][6] +S[2][20][2] )
                S[8][20][7] = gma0*S[2][20][7] + fak * (S[1][20][7] +S[2][7][7] )
                S[5][20][7] = gma1*S[2][20][7] + fak * (S[2][6][7] +S[2][20][4] )
                S[6][20][7] = gma2*S[2][20][7] + fak * (S[2][5][7] +S[2][20][3] )
                S[8][20][8] = gma0*S[2][20][8]  & 
                     +    fak * (S[1][20][8] +S[2][7][8] +2.d0*S[2][20][2] )
                S[5][20][8] = gma1*S[2][20][8] + fak * S[2][6][8]
                S[6][20][8] = gma2*S[2][20][8] + fak * S[2][5][8]
                S[8][20][9] = gma0*S[2][20][9] + fak * (S[1][20][9] +S[2][7][9] )
                S[5][20][9] = gma1*S[2][20][9] + fak * (S[2][6][9] +2.d0*S[2][20][3] )
                S[6][20][9] = gma2*S[2][20][9] + fak * S[2][5][9]
                S[8][20][10] = gma0*S[2][20][10] + fak * (S[1][20][10] +S[2][7][10] )
                S[5][20][10] = gma1*S[2][20][10] + fak * S[2][6][10]
                S[6][20][10] = gma2*S[2][20][10] + fak * (S[2][5][10] +2.d0*S[2][20][4])
                S[9][11][5] = gma1*S[3][11][5] + fak * (S[1][11][5] +S[3][11][2] )
                S[7][11][5] = gma2*S[3][11][5]
                S[9][11][6] = gma1*S[3][11][6] + fak * S[1][11][6]
                S[7][11][6] = gma2*S[3][11][6] + fak * S[3][11][2]
                S[9][11][7] = gma1*S[3][11][7] + fak * (S[1][11][7] +S[3][11][4] )
                S[7][11][7] = gma2*S[3][11][7] + fak * S[3][11][3]
                S[9][11][8] = gma1*S[3][11][8] + fak * S[1][11][8]
                S[7][11][8] = gma2*S[3][11][8]
                S[9][11][9] = gma1*S[3][11][9] + fak * (S[1][11][9] +2.d0*S[3][11][3] )
                S[7][11][9] = gma2*S[3][11][9]
                S[9][11][10] = gma1*S[3][11][10] + fak * S[1][11][10]
                S[7][11][10] = gma2*S[3][11][10] + fak2* S[3][11][4]
                S[9][12][5] = gma1*S[3][12][5]  & 
                     +    fak * (S[1][12][5] +3.d0*S[3][9][5] +S[3][12][2] )
                S[7][12][5] = gma2*S[3][12][5]
                S[9][12][6] = gma1*S[3][12][6] + fak * (S[1][12][6] +3.d0*S[3][9][6] )
                S[7][12][6] = gma2*S[3][12][6] + fak * S[3][12][2]
                S[9][12][7] = gma1*S[3][12][7]  & 
                     +    fak * (S[1][12][7] +3.d0*S[3][9][7] +S[3][12][4] )
                S[7][12][7] = gma2*S[3][12][7] + fak * S[3][12][3]
                S[9][12][8] = gma1*S[3][12][8] + fak * (S[1][12][8] +3.d0*S[3][9][8] )
                S[7][12][8] = gma2*S[3][12][8]
                S[9][12][9] = gma1*S[3][12][9]  & 
                     +    fak * (S[1][12][9] +3.d0*S[3][9][9] +2.d0*S[3][12][3] )
                S[7][12][9] = gma2*S[3][12][9]
                S[9][12][10] = gma1*S[3][12][10] + fak * (S[1][12][10]+3.d0*S[3][9][10])
                S[7][12][10] = gma2*S[3][12][10] + fak2* S[3][12][4]
                S[9][13][5] = gma1*S[3][13][5] + fak * (S[1][13][5] +S[3][13][2] )
                S[7][13][5] = gma2*S[3][13][5] + fak3* S[3][10][5]
                S[9][13][6] = gma1*S[3][13][6] + fak * S[1][13][6]
                S[7][13][6] = gma2*S[3][13][6] + fak * (3.d0*S[3][10][6] +S[3][13][2] )
                S[9][13][7] = gma1*S[3][13][7] + fak * (S[1][13][7] +S[3][13][4] )
                S[7][13][7] = gma2*S[3][13][7] + fak * (3.d0*S[3][10][7] +S[3][13][3] )
                S[9][13][8] = gma1*S[3][13][8] + fak * S[1][13][8]
                S[7][13][8] = gma2*S[3][13][8] + fak3* S[3][10][8]
                S[9][13][9] = gma1*S[3][13][9] + fak * (S[1][13][9] +2.d0*S[3][13][3] )
                S[7][13][9] = gma2*S[3][13][9] + fak3* S[3][10][9]
                S[9][13][10] = gma1*S[3][13][10] + fak * S[1][13][10]
                S[7][13][10] = gma2*S[3][13][10]  & 
                     +    fak * (3.d0*S[3][10][10] +2.d0*S[3][13][4] )
                S[9][14][5] = gma1*S[3][14][5] + fak *(S[1][14][5]+S[3][8][5]+S[3][14][2])
                S[7][14][5] = gma2*S[3][14][5]
                S[9][14][6] = gma1*S[3][14][6] + fak * (S[1][14][6] +S[3][8][6] )
                S[7][14][6] = gma2*S[3][14][6] + fak * S[3][14][2]
                S[9][14][7] = gma1*S[3][14][7] + fak *(S[1][14][7]+S[3][8][7]+S[3][14][4])
                S[7][14][7] = gma2*S[3][14][7] + fak * S[3][14][3]
                S[9][14][8] = gma1*S[3][14][8] + fak * (S[1][14][8] +S[3][8][8] )
                S[7][14][8] = gma2*S[3][14][8]
                S[9][14][9] = gma1*S[3][14][9]  & 
                     +    fak * (S[1][14][9] +S[3][8][9] +2.d0*S[3][14][3] )
                S[7][14][9] = gma2*S[3][14][9]
                S[9][14][10] = gma1*S[3][14][10] + fak * (S[1][14][10] +S[3][8][10] )
                S[7][14][10] = gma2*S[3][14][10] + fak2* S[3][14][4]
                S[9][15][5] = gma1*S[3][15][5]  & 
                     +    fak * (S[1][15][5] +2.d0*S[3][5][5] +S[3][15][2] )
                S[7][15][5] = gma2*S[3][15][5]
                S[9][15][6] = gma1*S[3][15][6] + fak * (S[1][15][6] +2.d0*S[3][5][6] )
                S[7][15][6] = gma2*S[3][15][6] + fak * S[3][15][2]
                S[9][15][7] = gma1*S[3][15][7]  & 
                     +    fak * (S[1][15][7] +2.d0*S[3][5][7] +S[3][15][4] )
                S[7][15][7] = gma2*S[3][15][7] + fak * S[3][15][3]
                S[9][15][8] = gma1*S[3][15][8] + fak * (S[1][15][8] +2.d0*S[3][5][8] )
                S[7][15][8] = gma2*S[3][15][8]
                S[9][15][9] = gma1*S[3][15][9]  & 
                     +    fak * (S[1][15][9] +2.d0*S[3][5][9] +2.d0*S[3][15][3] )
                S[7][15][9] = gma2*S[3][15][9]
                S[9][15][10] = gma1*S[3][15][10] + fak * (S[1][15][10]+2.d0*S[3][5][10])
                S[7][15][10] = gma2*S[3][15][10] + fak2* S[3][15][4]
                S[9][16][5] = gma1*S[3][16][5] + fak * (S[1][16][5] +S[3][16][2] )
                S[7][16][5] = gma2*S[3][16][5] + fak * S[3][8][5]
                S[9][16][6] = gma1*S[3][16][6] + fak * S[1][16][6]
                S[7][16][6] = gma2*S[3][16][6] + fak * (S[3][8][6] +S[3][16][2] )
                S[9][16][7] = gma1*S[3][16][7] + fak * (S[1][16][7] +S[3][16][4] )
                S[7][16][7] = gma2*S[3][16][7] + fak * (S[3][8][7] +S[3][16][3] )
                S[9][16][8] = gma1*S[3][16][8] + fak * S[1][16][8]
                S[7][16][8] = gma2*S[3][16][8] + fak * S[3][8][8]
                S[9][16][9] = gma1*S[3][16][9] + fak * (S[1][16][9] +2.d0*S[3][16][3] )
                S[7][16][9] = gma2*S[3][16][9] + fak * S[3][8][9]
                S[9][16][10] = gma1*S[3][16][10] + fak * S[1][16][10]
                S[7][16][10] = gma2*S[3][16][10] + fak * (S[3][8][10] +2.d0*S[3][16][4])
                S[9][17][5] = gma1*S[3][17][5] + fak * (S[1][17][5] +S[3][17][2] )
                S[7][17][5] = gma2*S[3][17][5] + fak2* S[3][6][5]
                S[9][17][6] = gma1*S[3][17][6] + fak * S[1][17][6]
                S[7][17][6] = gma2*S[3][17][6] + fak * (2.d0*S[3][6][6] +S[3][17][2] )
                S[9][17][7] = gma1*S[3][17][7] + fak * (S[1][17][7] +S[3][17][4] )
                S[7][17][7] = gma2*S[3][17][7] + fak * (2.d0*S[3][6][7] +S[3][17][3] )
                S[9][17][8] = gma1*S[3][17][8] + fak * S[1][17][8]
                S[7][17][8] = gma2*S[3][17][8] + fak2* S[3][6][8]
                S[9][17][9] = gma1*S[3][17][9] + fak * (S[1][17][9] +2.d0*S[3][17][3] )
                S[7][17][9] = gma2*S[3][17][9] + fak2* S[3][6][9]
                S[9][17][10] = gma1*S[3][17][10] + fak * S[1][17][10]
                S[7][17][10] = gma2*S[3][17][10] + fak2* (S[3][6][10] +S[3][17][4] )
                S[9][18][5] = gma1*S[3][18][5]  & 
                     +    fak * (S[1][18][5] +2.d0*S[3][7][5] +S[3][18][2] )
                S[7][18][5] = gma2*S[3][18][5] + fak * S[3][9][5]
                S[9][18][6] = gma1*S[3][18][6] + fak * (S[1][18][6] +2.d0*S[3][7][6] )
                S[7][18][6] = gma2*S[3][18][6] + fak * (S[3][9][6] +S[3][18][2] )
                S[9][18][7] = gma1*S[3][18][7]  & 
                     +    fak * (S[1][18][7] +2.d0*S[3][7][7] +S[3][18][4] )
                S[7][18][7] = gma2*S[3][18][7] + fak * (S[3][9][7] +S[3][18][3] )
                S[9][18][8] = gma1*S[3][18][8] + fak * (S[1][18][8] +2.d0*S[3][7][8] )
                S[7][18][8] = gma2*S[3][18][8] + fak * S[3][9][8]
                S[9][18][9] = gma1*S[3][18][9]  & 
                     +    fak * (S[1][18][9] +2.d0*S[3][7][9] +2.d0*S[3][18][3] )
                S[7][18][9] = gma2*S[3][18][9] + fak * S[3][9][9]
                S[9][18][10] = gma1*S[3][18][10] + fak * (S[1][18][10]+2.d0*S[3][7][10])
                S[7][18][10] = gma2*S[3][18][10] + fak * (S[3][9][10] +2.d0*S[3][18][4])
                S[9][19][5] = gma1*S[3][19][5] + fak*(S[1][19][5]+S[3][10][5]+S[3][19][2])
                S[7][19][5] = gma2*S[3][19][5] + fak2* S[3][7][5]
                S[9][19][6] = gma1*S[3][19][6] + fak * (S[1][19][6] +S[3][10][6] )
                S[7][19][6] = gma2*S[3][19][6] + fak * (2.d0*S[3][7][6] +S[3][19][2] )
                S[9][19][7] = gma1*S[3][19][7] + fak*(S[1][19][7]+S[3][10][7]+S[3][19][4])
                S[7][19][7] = gma2*S[3][19][7] + fak * (2.d0*S[3][7][7] +S[3][19][3] )
                S[9][19][8] = gma1*S[3][19][8] + fak * (S[1][19][8] +S[3][10][8] )
                S[7][19][8] = gma2*S[3][19][8] + fak2* S[3][7][8]
                S[9][19][9] = gma1*S[3][19][9]  & 
                     +    fak * (S[1][19][9] +S[3][10][9] +2.d0*S[3][19][3] )
                S[7][19][9] = gma2*S[3][19][9] + fak2* S[3][7][9]
                S[9][19][10] = gma1*S[3][19][10] + fak * (S[1][19][10] +S[3][10][10] )
                S[7][19][10] = gma2*S[3][19][10] + fak2* (S[3][7][10] +S[3][19][4] )
                S[9][20][5] = gma1*S[3][20][5]  & 
                     +    fak * (S[1][20][5] +S[3][6][5] + S[3][20][2] )
                S[7][20][5] = gma2*S[3][20][5] + fak * S[3][5][5]
                S[9][20][6] = gma1*S[3][20][6] + fak * (S[1][20][6] +S[3][6][6] )
                S[7][20][6] = gma2*S[3][20][6] + fak * (S[3][5][6] +S[3][20][2] )
                S[9][20][7] = gma1*S[3][20][7] + fak *(S[1][20][7]+S[3][6][7]+S[3][20][4])
                S[7][20][7] = gma2*S[3][20][7] + fak * (S[3][5][7] +S[3][20][3] )
                S[9][20][8] = gma1*S[3][20][8] + fak * (S[1][20][8] +S[3][6][8] )
                S[7][20][8] = gma2*S[3][20][8] + fak * S[3][5][8]
                S[9][20][9] = gma1*S[3][20][9]  & 
                     +    fak * (S[1][20][9] +S[3][6][9] +2.d0*S[3][20][3] )
                S[7][20][9] = gma2*S[3][20][9] + fak * S[3][5][9]
                S[9][20][10] = gma1*S[3][20][10] + fak * (S[1][20][10] +S[3][6][10] )
                S[7][20][10] = gma2*S[3][20][10] + fak * (S[3][5][10] +2.d0*S[3][20][4])
                S[10][11][5] = gma2*S[4][11][5] + fak * S[1][11][5]
                S[10][11][6] = gma2*S[4][11][6] + fak * (S[1][11][6] +S[4][11][2] )
                S[10][11][7] = gma2*S[4][11][7] + fak * (S[1][11][7] +S[4][11][3] )
                S[10][11][8] = gma2*S[4][11][8] + fak * S[1][11][8]
                S[10][11][9] = gma2*S[4][11][9] + fak * S[1][11][9]
                S[10][11][10] = gma2*S[4][11][10] + fak *(S[1][11][10]+2.d0*S[4][11][4])
                S[10][12][5] = gma2*S[4][12][5] + fak * S[1][12][5]
                S[10][12][6] = gma2*S[4][12][6] + fak * (S[1][12][6] +S[4][12][2] )
                S[10][12][7] = gma2*S[4][12][7] + fak * (S[1][12][7] +S[4][12][3] )
                S[10][12][8] = gma2*S[4][12][8] + fak * S[1][12][8]
                S[10][12][9] = gma2*S[4][12][9] + fak * S[1][12][9]
                S[10][12][10] = gma2*S[4][12][10] + fak *(S[1][12][10]+2.d0*S[4][12][4])
                S[10][13][5] = gma2*S[4][13][5] + fak * (S[1][13][5] +3.d0*S[4][10][5] )
                S[10][13][6] = gma2*S[4][13][6]  & 
                     +    fak * (S[1][13][6] +3.d0*S[4][10][6] +S[4][13][2] )
                S[10][13][7] = gma2*S[4][13][7]  & 
                     +    fak * (S[1][13][7] +3.d0*S[4][10][7] +S[4][13][3] )
                S[10][13][8] = gma2*S[4][13][8] + fak * (S[1][13][8] +3.d0*S[4][10][8] )
                S[10][13][9] = gma2*S[4][13][9] + fak * (S[1][13][9] +3.d0*S[4][10][9] )
                S[10][13][10] = gma2*S[4][13][10]  & 
                     +    fak * (S[1][13][10] +3.d0*S[4][10][10] +2.d0*S[4][13][4] )
                S[10][14][5] = gma2*S[4][14][5] + fak * S[1][14][5]
                S[10][14][6] = gma2*S[4][14][6] + fak * (S[1][14][6] +S[4][14][2] )
                S[10][14][7] = gma2*S[4][14][7] + fak * (S[1][14][7] +S[4][14][3] )
                S[10][14][8] = gma2*S[4][14][8] + fak * S[1][14][8]
                S[10][14][9] = gma2*S[4][14][9] + fak * S[1][14][9]
                S[10][14][10] = gma2*S[4][14][10] + fak *(S[1][14][10]+2.d0*S[4][14][4])
                S[10][15][5] = gma2*S[4][15][5] + fak * S[1][15][5]
                S[10][15][6] = gma2*S[4][15][6] + fak * (S[1][15][6] +S[4][15][2] )
                S[10][15][7] = gma2*S[4][15][7] + fak * (S[1][15][7] +S[4][15][3] )
                S[10][15][8] = gma2*S[4][15][8] + fak * S[1][15][8]
                S[10][15][9] = gma2*S[4][15][9] + fak * S[1][15][9]
                S[10][15][10] = gma2*S[4][15][10] + fak *(S[1][15][10]+2.d0*S[4][15][4])
                S[10][16][5] = gma2*S[4][16][5] + fak * (S[1][16][5] +S[4][8][5] )
                S[10][16][6] = gma2*S[4][16][6] + fak*(S[1][16][6]+S[4][8][6]+S[4][16][2])
                S[10][16][7] = gma2*S[4][16][7] + fak*(S[1][16][7]+S[4][8][7]+S[4][16][3])
                S[10][16][8] = gma2*S[4][16][8] + fak * (S[1][16][8] +S[4][8][8] )
                S[10][16][9] = gma2*S[4][16][9] + fak * (S[1][16][9] +S[4][8][9] )
                S[10][16][10] = gma2*S[4][16][10]  & 
                     +    fak * (S[1][16][10] +S[4][8][10] +2.d0*S[4][16][4] )
                S[10][17][5] = gma2*S[4][17][5] + fak * (S[1][17][5] +2.d0*S[4][6][5] )
                S[10][17][6] = gma2*S[4][17][6]  & 
                     +    fak * (S[1][17][6] +2.d0*S[4][6][6] +S[4][17][2] )
                S[10][17][7] = gma2*S[4][17][7]  & 
                     +    fak * (S[1][17][7] +2.d0*S[4][6][7] +S[4][17][3] )
                S[10][17][8] = gma2*S[4][17][8] + fak * (S[1][17][8] +2.d0*S[4][6][8] )
                S[10][17][9] = gma2*S[4][17][9] + fak * (S[1][17][9] +2.d0*S[4][6][9] )
                S[10][17][10] = gma2*S[4][17][10]  & 
                     +    fak * (S[1][17][10] +2.d0*S[4][6][10] +2.d0*S[4][17][4] )
                S[10][18][5] = gma2*S[4][18][5] + fak * (S[1][18][5] +S[4][9][5] )
                S[10][18][6] = gma2*S[4][18][6] + fak*(S[1][18][6]+S[4][9][6]+S[4][18][2])
                S[10][18][7] = gma2*S[4][18][7] + fak*(S[1][18][7]+S[4][9][7]+S[4][18][3])
                S[10][18][8] = gma2*S[4][18][8] + fak * (S[1][18][8] +S[4][9][8] )
                S[10][18][9] = gma2*S[4][18][9] + fak * (S[1][18][9] +S[4][9][9] )
                S[10][18][10] = gma2*S[4][18][10]  & 
                     +    fak * (S[1][18][10] +S[4][9][10] +2.d0*S[4][18][4] )
                S[10][19][5] = gma2*S[4][19][5] + fak * (S[1][19][5] +2.d0*S[4][7][5] )
                S[10][19][6] = gma2*S[4][19][6]  & 
                     +    fak * (S[1][19][6] +2.d0*S[4][7][6] +S[4][19][2] )
                S[10][19][7] = gma2*S[4][19][7]  & 
                     +    fak * (S[1][19][7] +2.d0*S[4][7][7] +S[4][19][3] )
                S[10][19][8] = gma2*S[4][19][8] + fak * (S[1][19][8] +2.d0*S[4][7][8] )
                S[10][19][9] = gma2*S[4][19][9] + fak * (S[1][19][9] +2.d0*S[4][7][9] )
                S[10][19][10] = gma2*S[4][19][10]  & 
                     +    fak * (S[1][19][10] +2.d0*S[4][7][10] +2.d0*S[4][19][4] )
                S[10][20][5] = gma2*S[4][20][5] + fak * (S[1][20][5] +S[4][5][5] )
                S[10][20][6] = gma2*S[4][20][6] + fak*(S[1][20][6]+S[4][5][6]+S[4][20][2])
                S[10][20][7] = gma2*S[4][20][7] + fak*(S[1][20][7]+S[4][5][7]+S[4][20][3])
                S[10][20][8] = gma2*S[4][20][8] + fak * (S[1][20][8] +S[4][5][8] )
                S[10][20][9] = gma2*S[4][20][9] + fak * (S[1][20][9] +S[4][5][9] )
                S[10][20][10] = gma2*S[4][20][10]  & 
                     +    fak * (S[1][20][10] +S[4][5][10] +2.d0*S[4][20][4] )
             } */

            // data is now stored in unnormalized cartesian Gaussians in the multiarray
            // Now, weird-looking construction since multiarray is not accessible for ub::prod
            //              s  px  py  pz dxz dyz dxy d3z2-r2 dx2-y2  f1  f2  f3  f4  f5  f6  f7
            int istart[] = {0, 1, 2, 3, 5, 6, 4, 7, 7, 12, 10, 11, 11, 10, 19, 15};
            int istop[]  = {0, 1, 2, 3, 5, 6, 4, 9, 8, 17, 16, 18, 13, 14, 19, 17};

            // ub::vector<ub::matrix<double> >& _subvector
            // which ones do we want to store
            int _offset_gw = _shell_gw->getOffset();
            int _offset_alpha = _shell_alpha->getOffset();
            int _offset_gamma = _shell_gamma->getOffset();

            // prepare transformation matrices
            int _ntrafo_gw = _shell_gw->getNumFunc() + _offset_gw;
            int _ntrafo_alpha = _shell_alpha->getNumFunc() + _offset_alpha;
            int _ntrafo_gamma = _shell_gamma->getNumFunc() + _offset_gamma;

            ub::matrix<double> _trafo_gw = ub::zero_matrix<double>(_ntrafo_gw, _ngw);
            ub::matrix<double> _trafo_alpha = ub::zero_matrix<double>(_ntrafo_alpha, _nalpha);
            ub::matrix<double> _trafo_gamma = ub::zero_matrix<double>(_ntrafo_gamma, _ngamma);

            // get transformation matrices
            this->getTrafo(_trafo_gw, _lmax_gw, _decay_gw);
            this->getTrafo(_trafo_alpha, _lmax_alpha, _decay_alpha);
            this->getTrafo(_trafo_gamma, _lmax_gamma, _decay_gamma);

            // transform from unnormalized cartesians to normalized sphericals
            // container with indices starting at zero
            ma_type S_sph;
            S_sph.resize(extents[ _ntrafo_alpha ][ _ntrafo_gw  ][  _ntrafo_gamma ]);

            for (int _i_gw = 0; _i_gw < _ntrafo_gw; _i_gw++) {
                for (int _i_alpha = 0; _i_alpha < _ntrafo_alpha; _i_alpha++) {
                    for (int _i_gamma = 0; _i_gamma < _ntrafo_gamma; _i_gamma++) {
             
                        S_sph[ _i_alpha ][ _i_gw ][ _i_gamma ] = 0.0;
                        
                        for (int _i_gw_t = istart[ _i_gw ]; _i_gw_t <= istop[ _i_gw ]; _i_gw_t++) {
                            for (int _i_alpha_t = istart[ _i_alpha ]; _i_alpha_t <= istop[ _i_alpha ]; _i_alpha_t++) {
                                for (int _i_gamma_t = istart[ _i_gamma ]; _i_gamma_t <= istop[ _i_gamma ]; _i_gamma_t++) {

                                     S_sph[ _i_alpha ][ _i_gw ][ _i_gamma ] += S[ _i_alpha_t +1 ][ _i_gw_t +1 ][ _i_gamma_t +1]  
                                                * _trafo_alpha(_i_alpha, _i_alpha_t) * _trafo_gw(_i_gw, _i_gw_t) * _trafo_gamma(_i_gamma, _i_gamma_t);
                                    

                                }
                            }
                        }
                    }
                }
            }

            // only store the parts, we need
            for (int _i_gw = 0; _i_gw < _shell_gw->getNumFunc(); _i_gw++) {
                for (int _i_alpha = 0; _i_alpha < _shell_alpha->getNumFunc(); _i_alpha++) {
                    for (int _i_gamma = 0; _i_gamma < _shell_gamma->getNumFunc(); _i_gamma++) {

                
                         int _i_index = _shell_gamma->getNumFunc() *  _i_gw   + _i_gamma;
                         
                         _subvector( _i_alpha, _i_index ) = S_sph[ _offset_alpha + _i_alpha ][ _offset_gw + _i_gw ][ _offset_gamma + _i_gamma ];
                        
                    }
                }
            }
            
 
            return true;


        }
        
        
        
        
        void TCMatrix::Print(string _ident) {
            cout << "\n" << endl;
            for (int k = 0; k < this->mtotal; k++){
                    for (unsigned int i = 0; i< _matrix[1].size1() ; i++) {
                        for (int j = 0; j< this->ntotal; j++) {
                           cout << _ident << "[" << i+1 << ":" << k + 1 << ":" << j + 1 << "] " << this->_matrix[k](i, j) << endl;
                        }
                }
            }
        }

        void TCMatrix::getTrafo(ub::matrix<double>& _trafo, int _lmax, const double& _decay) {
            // s-functions
            _trafo(0, 0) = 1.0; // s

            // p-functions
            if (_lmax > 0) {
                //cout << _trafo_row.size1() << ":" << _trafo_row.size2() << endl;
                _trafo(1, 1) = 2.0 * sqrt(_decay);
                _trafo(2, 2) = 2.0 * sqrt(_decay);
                _trafo(3, 3) = 2.0 * sqrt(_decay);
            }

            // d-functions
            if (_lmax > 1) {
                _trafo(4, 5) = 4.0 * _decay; // dxz
                _trafo(5, 6) = _trafo(4, 5); // dyz
                _trafo(6, 4) = _trafo(4, 5); // dxy
                _trafo(7, 7) = -2.0 * _decay / sqrt(3.0); // d3z2-r2 (dxx)
                _trafo(7, 8) = _trafo(7, 7); // d3z2-r2 (dyy)
                _trafo(7, 9) = -2.0 * _trafo(7, 7); // d3z2-r2 (dzz)
                _trafo(8, 7) = 2.0 * _decay; // dx2-y2 (dxx)
                _trafo(8, 8) = -_trafo(8, 7); // dx2-y2 (dzz)
            }

            // f-functions
            if (_lmax > 2) {
                _trafo(9, 12) = 4.0 * 2.0 * pow(_decay, 1.5); // f1 (f??)
                _trafo(9, 15) = -1.5 * _trafo(9, 12); // f1 (f??)
                _trafo(9, 17) = _trafo(9, 15); // f1 (f??)

                _trafo(10, 16) = 4.0 * 2.0 * sqrt(2.0) / sqrt(5.0) * pow(_decay, 1.5); // f2 (f??)
                _trafo(10, 10) = -0.25 * _trafo(10, 16); // f2 f(??)
                _trafo(10, 14) = _trafo(10, 10); // f2 f(??)

                _trafo(11, 18) = _trafo(10, 16); // f3 (f??)
                _trafo(11, 13) = -0.25 * _trafo(11, 18); // f3 f(??)
                _trafo(11, 11) = _trafo(11, 13); // f3 f(??)            

                _trafo(12, 13) = 3.0 * 2.0 * sqrt(2.0) / sqrt(3.0) * pow(_decay, 1.5); // f4 (f??)
                _trafo(12, 11) = -_trafo(12, 13) / 3.0; // f4 (f??)

                _trafo(13, 10) = -_trafo(12, 11); // f5 (f??)
                _trafo(13, 14) = -_trafo(12, 13); // f5 (f??)

                _trafo(14, 19) = 8.0 * pow(_decay, 1.5); // f6 (f??)

                _trafo(15, 15) = 0.5 * _trafo(14, 19); // f7 (f??)
                _trafo(15, 17) = -_trafo(15, 15); // f7 (f??)
            }

            // g-functions
            if (_lmax > 3) {
                _trafo(16, 22) = 8.0 * 2.0 / sqrt(105.0) * pow(_decay, 2.0);
                _trafo(16, 21) = 3.0 * 2.0 / sqrt(105.0) * pow(_decay, 2.0);
                _trafo(16, 20) = _trafo(16, 21);
                _trafo(16, 29) = -3.0 * _trafo(16, 22);
                _trafo(16, 31) = 2.0 * _trafo(16, 21);
                _trafo(16, 30) = _trafo(16, 29);
                _trafo(16, 5) = _trafo(16, 31);

                /* vv(17,:) =  (/   23,  22, 21, 30, 32, 31,   6 /) ! g
                   cc(17,:) =  (/    8,  3, 3, -24, 6, -24,    6 /)
                   normConst(17,:) = (/ 2.d0/sqrt(105.d0) ,2.d0  /)
                 */
                _trafo(17, 26) = 4.0 * 4.0 * sqrt(2.0) / sqrt(21.0) * pow(_decay, 2.0);
                _trafo(17, 25) = -0.75 * _trafo(17, 26);
                _trafo(17, 33) = _trafo(17, 25);

                /* vv(18,:) =  (/   27,  26, 34,  0,  0,  0,   3 /) ! g
                   cc(18,:) =  (/    4,  -3, -3,  0,  0,  0,   3 /)
                   normConst(18,:) = (/ 4.d0*sqrt(2.d0)/sqrt(21.d0) ,2.d0  /)
                 */

                _trafo(18, 28) = _trafo(17, 26);
                _trafo(18, 32) = _trafo(17, 25);
                _trafo(18, 27) = _trafo(17, 25);

                /* vv(19,:) =  (/   29,  33, 28,  0,  0,  0,   3 /) ! g 
                   cc(19,:) =  (/    4,  -3, -3,  0,  0,  0,   3 /)
                   normConst(19,:) = (/ 4.d0*sqrt(2.d0)/sqrt(21.d0) ,2.d0  /)
                 */

                _trafo(19, 34) = 6.0 * 8.0 / sqrt(21.0) * pow(_decay, 2.0);
                _trafo(19, 23) = -_trafo(19, 34) / 6.0;
                _trafo(19, 24) = _trafo(19, 23);

                /* vv(20,:) =  (/   35,  24, 25,  0,  0,  0,   3 /) ! g
                   cc(20,:) =  (/    6,  -1, -1,  0,  0,  0,   3 /)
                   normConst(20,:) = (/ 8.d0/sqrt(21.d0) ,2.d0  /)
                 */

                _trafo(20, 29) = 6.0 * 4.0 / sqrt(21.0) * pow(_decay, 2.0);
                _trafo(20, 20) = -_trafo(20, 29) / 6.0;
                _trafo(20, 30) = -_trafo(20, 29);
                _trafo(20, 21) = -_trafo(20, 20);

                /* vv(21,:) =  (/   30,  21, 31, 22,  0,  0,   4 /) ! g
                   cc(21,:) =  (/    6,  -1, -6, 1,  0,  0,    4 /)
                   normConst(21,:) = (/ 4.d0/sqrt(21.d0) ,2.d0  /)
                 */

                _trafo(21, 25) = 4.0 * sqrt(2.0) / sqrt(3.0) * pow(_decay, 2.0);
                _trafo(21, 33) = -3.0 * _trafo(21, 25);

                /* vv(22,:) =  (/   26,  34,  0,  0,  0,  0,   2 /) ! g
                   cc(22,:) =  (/    1,  -3,  0,  0,  0,  0,   2 /)
                   normConst(22,:) = (/ 4.d0*sqrt(2.d0)/sqrt(3.d0) ,2.d0  /)
                 */

                _trafo(22, 32) = -_trafo(21, 33);
                _trafo(22, 27) = -_trafo(21, 25);

                /* vv(23,:) =  (/   33,  28,  0,  0,  0,  0,   2 /) ! g
                   cc(23,:) =  (/    3,  -1,  0,  0,  0,  0,   2 /)
                   normConst(23,:) = (/ 4.d0*sqrt(2.d0)/sqrt(3.d0) ,2.d0  /)
                 */

                _trafo(23, 23) = 8.0 / sqrt(3.0) * pow(_decay, 2.0);
                _trafo(23, 24) = -_trafo(23, 23);

                /* vv(24,:) =  (/   24,  25,  0,  0,  0,  0,   2 /) ! g 
                   cc(24,:) =  (/    1,  -1,  0,  0,  0,  0,   2 /)
                   normConst(24,:) = (/ 8.d0/sqrt(3.d0) ,2.d0  /)
                 */

                _trafo(24, 20) = 2.0 / sqrt(3.0) * pow(_decay, 2.0);
                _trafo(24, 21) = _trafo(24, 20);
                _trafo(24, 31) = -6.0 * _trafo(24, 20);

                /* vv(25,:) =  (/   21,  22, 32,  0,  0,  0,   3 /) ! g
                   cc(25,:) =  (/    1,  1, -6,  0,  0,  0,   3  /)
                   normConst(25,:) = (/ 2.d0/sqrt(3.d0) ,2.d0  /)
                 */

            }


        }

        int TCMatrix::getBlockSize(int _lmax) {
            int _block_size;
            if (_lmax == 0) {
                _block_size = 1;
            } // s
            if (_lmax == 1) {
                _block_size = 4;
            } // p
            if (_lmax == 2) {
                _block_size = 10;
            } // d
            if (_lmax == 3) {
                _block_size = 20;
            } // f
            if (_lmax == 4) {
                _block_size = 35;
            } // g

            return _block_size;
        }



    }
}

