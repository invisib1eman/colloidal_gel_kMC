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
    double BoxLength;
    double R_hardcore;
    double well_depth = -10000;
    double cm_L;
    double R_hardcore_DH;
    double cutoff_distance;
    double morse_well_depth;
    double morse_a;
    double morse_r0;
    double yukawa_e;
    double yukawa_debye_length;
    double morse_potential(double r2)
    {
        double r = sqrt(r2);
        double shift_morse = morse_well_depth*(exp(-2*(cutoff_distance-morse_r0)/morse_a)-2*exp(-(cutoff_distance-morse_r0)/morse_a));
        if (r > cutoff_distance)
        {
            return 0;
        }
        else
        {
            double U = morse_well_depth*(exp(-2*(r-morse_r0)*morse_a)-2*exp(-(r-morse_r0)*morse_a)) - shift_morse;
            return U;
        }
    }
    double Debye_Huckel(double r2)
    {
        
        double r = sqrt(r2);
        if (r > cutoff_distance)
        {
            return 0;
        }
        else
        {
            double prefactor = charge*charge*bjerrum_length*exp(2*R_hardcore_DH/debye_length)/((1+R_hardcore_DH/debye_length)*(1+R_hardcore_DH/debye_length));
            double shift = prefactor*exp(-cutoff_distance/debye_length)/cutoff_distance;
            double U = prefactor*exp(-r/debye_length)/r - shift;
            return U;
        }
    }
    double yukawa(double r2)
    {
        double r = sqrt(r2);
        if (r > cutoff_distance)
        {
            return 0;
        }
        else
        {
            return yukawa_e*cm_L*exp(-(r-cm_L)/yukawa_debye_length)/r;
        }
    }
    double potential_well(double r2)
    {
        double r = sqrt(r2);
        if (r < cm_L && r > 2*R_hardcore)
        {
            return well_depth;
        }
        else if (r <= 2*R_hardcore)
        {
            return 100000;
        }
        else
        {
            return 0;
        }
    }
    double total_energy_yukawa(double r2)
    {
        return yukawa(r2)+potential_well(r2);
    }
    double total_energy(double r2)
    {
        return Debye_Huckel(r2)+potential_well(r2);
    }
};
#endif

