//MIXING: enegies.h Energy Class (Revision Date: Feb 2, 2012)
#ifndef _ENERGY_H
#define _ENERGY_H
#include "system.h"
#include "utils.h"
#include "particle.h"
#include "xyz.h"
#include "particle.h"
class Energy
{
public:
    double debye_length;
    double bjerrum_length;
    double charge;
    double L;
    double R_hardcore;
    double well_depth = -10000;
    double well_width;
    double well_edge;
    double Debye_Huckel(double r2)
    {
        double r = sqrt(r2);
        double prefactor = charge*charge*bjerrum_length*exp(2*R_hardcore/debye_length)/((1+R_hardcore/debye_length)*(1+R_hardcore/debye_length));
        double U = prefactor*exp(-r/debye_length)/r;
        return U;
    }
    double potential_well(double r2)
    {
        double r = sqrt(r2);
        if (r < 2*well_edge && r > 2*R_hardcore)
        {
            return well_depth;
        }
        else if (r <= 2*R_hardcore)
        {
            return 10000;
        }
        else
        {
            return 0;
        }
    }
    double total_energy(double r2)
    {
        return Debye_Huckel(r2)+potential_well(r2);
    }
};
#endif

