#ifndef PARAMETERS_H
#define PARAMETERS_H

struct bag_site{

	int bag_num;
	bool n0;
	bool n1;
	bool n2;
	bool n3;
};
 
namespace param
{
    // forward declarations only
    extern const int Nt;
    extern const int Ns;
		extern const int V;
		extern double m;
		extern const int nmeas;
		extern const int nequi;
		extern const int nskip;
		extern int max_label;
}
 
#endif
