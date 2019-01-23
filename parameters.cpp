#include "parameters.h"
#include<chrono>

namespace param
{
    //Simulation parameters
    extern const int Nt(4);					//lattice size in temporal direction
    extern const int Ns(4);					//lattice siye in spatial direction
		extern const int V(Nt*Ns);			//# of lattice sites, i.e. the volume
		extern double m(1.0);			//bare mass
		extern const int nmeas(1e4);		//m
		extern const int nequi(1e4);
		extern const int nskip(1e2);
		extern int max_label(0);
}
