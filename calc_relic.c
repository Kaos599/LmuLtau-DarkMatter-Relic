#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "micromegas.h"
#include "micromegas_aux.h"
#include "libg.h"

/*
 * C wrapper for micrOMEGAs.
 * Takes mchi, mzpl, gp as command line arguments, calculates relic density Omega_h2 
 * and thermally averaged cross-section <sigma v> today.
 */
int main(int argc, char** argv) {
    int err;
    double Xf;
    double Omega;
    
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <mchi> <mzpl> <gp>\n", argv[0]);
        return 1;
    }
    
    double mchi = atof(argv[1]);
    double mzpl = atof(argv[2]);
    double gp   = atof(argv[3]);
    
    // Assign parameters to micrOMEGAs engine
    assignValW("Mchi", mchi);
    assignValW("MZp", mzpl);
    assignValW("gp", gp);
    
    // Sort odd particles to find the dark matter candidate
    err = sortOddParticles("~chi");
    if (err) {
        fprintf(stderr, "Error in sortOddParticles: %d\n", err);
        return 1;
    }
    
    // Calculate relic density
    int err_omega = 0;
    Omega = darkOmega(&Xf, 1, 1e-5, &err_omega);
    if (err_omega) {
        fprintf(stderr, "Error in darkOmega: %d\n", err_omega);
        return 1;
    }
    
    // Output parsed values for Python automation script
    printf("OMEGA_H2=%.8e\n", Omega);
    printf("XF=%.8e\n", Xf);
    
    // Thermally averaged cross section at freeze out: Tf = Mchi / Xf
    double Tf = mchi / Xf;
    double sigmav_fo = vSigma(Tf, 1e-5, &err);
    printf("SIGMAV_FO=%.8e\n", sigmav_fo);
    
    // Thermally averaged cross section today (v ~ 10^-3)
    // T_today ~ Mchi * (v0^2 / 6) ~ Mchi * 1.67e-7 (sub-GeV to TeV scale)
    double T_today = mchi * 1.67e-7;
    double sigmav_today = vSigma(T_today, 1e-5, &err);
    printf("SIGMAV_TODAY=%.8e\n", sigmav_today);
    
    return 0;
}
