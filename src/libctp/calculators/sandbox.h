#ifndef SANDBOX_H
#define SANDBOX_H

#include <votca/ctp/qmcalculator.h>

namespace votca { namespace ctp {

class Sandbox : public QMCalculator
{
public:
    Sandbox() { };
   ~Sandbox() { };

    void Initialize(QMTopology *top, Property *opt);
    bool EvaluateFrame(QMTopology *top);

private:
    int     _ID;

};

}} /* exit namespace votca::ctp */

#endif /* SANDBOX_H */
