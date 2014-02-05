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

#ifndef VOTCA_CTP_TRAJ_IO_H
#define VOTCA_CTP_TRAJ_IO_H

#include <iostream>
#include <votca/ctp/topology.h>

namespace votca { namespace ctp {

class trajIO
{
public:
    trajIO(){}
    virtual ~trajIO(){}
    virtual void read (std::string,Topology*) = 0;
    virtual void write() = 0;
private:
    
};

} /*namespace votca END */ } /* namespace ctp END */
#endif /* VOTCA_CTP_TRAJ_IO_H */