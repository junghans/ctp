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


#include <votca/ctp/extractorfactory.h>
#include "votca_config.h"

#include "extractors/energyextractor.h"
#include "extractors/integralsextractor.h"
#include "extractors/ratesextractor.h"
#include "extractors/trajextractor.h"
#include "extractors/segmentsextractor.h"
#include "extractors/pairsextractor.h"
#include "extractors/occupationsextractor.h"



namespace votca { namespace ctp {

void ExtractorFactory::RegisterAll(void)
{	
        Extractors().Register<EnergyExtractor>             ("energy");
        Extractors().Register<IntegralsExtractor>          ("integrals");
        Extractors().Register<RatesExtractor>              ("rates");
        Extractors().Register<OccupationsExtractor>        ("occupations");
        Extractors().Register<TrajExtractor>               ("trajectory");
        Extractors().Register<SegmentsExtractor>           ("segments");
        Extractors().Register<PairsExtractor>              ("pairs");
}

}}