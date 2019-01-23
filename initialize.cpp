#include "initialize.h"

const unsigned seed1(std::chrono::system_clock::now().time_since_epoch().count());
//std::default_random_engine generator(seed);
std::ranlux48 generator1(seed1);

using namespace std;


void neib_init (const int &N1, const int &N2, int neib[][4])
{
    int k,xm,xp,yp,ym;

    for(int x=0; x < N1; x++)
    {
        xp= x+1;
        xm= x-1;

        if (xp == N1) xp = 0;
        if (xm == -1) xm = N1 - 1;

        for(int y=0; y < N2; y++)
        {
            yp = y + 1;
            ym = y - 1;

            if (yp == N2) yp = 0;
            if (ym == -1) ym = N2 - 1;

            k = x + y*N1;

            neib[k][0]= x + yp*N1;
            neib[k][1]= xp + y*N1;
            neib[k][2]= x + ym*N1;
            neib[k][3]= xm + y*N1;
        }
    }

    return;
}

//this function initalizes a monomer dimer system containing only (2,2,1) configs, only on even lattices
void init_md_config(const int neib[][4], vector<int>& site_type, vector<vector<int>>& dimer, vector<int>& bag_vec)
{
	for(int i = 0; i < param::V; i++)
	{
		bag_vec.at(i) = 0; //no bags
		dimer.at(i).at(0) = 0; //only dimers in spatial direction
		if(site_type.at(i) < 0 && site_type.at(neib[i][1]) < 0)
		{
			site_type.at(i) = 2;
			site_type.at(neib[i][1]) = 2;
			dimer.at(i).at(1) = 1;
		}
		else dimer.at(i).at(0) = 0;
	}
}

//initialize the lattice as a single bag
void config_init(vector<int>& site_type, vector<vector<int>>& dimer, vector<int>& bag_vec, vector<vector<bool>>& bag_site_vec)
{
	for(int i = 0; i < param::V; i++)
	{
		site_type.at(i) = 5; //site type 5 = bag!
		bag_vec.at(i) = 1;
		bag_site_vec.at(i).at(0) = true;
		bag_site_vec.at(i).at(1) = true;
		bag_site_vec.at(i).at(2) = true;
		bag_site_vec.at(i).at(3) = true;
		dimer.at(i).at(0) = 0;
		dimer.at(i).at(1) = 0;
	}
}


//void sig_init (vector<vector<int>>& sig_link)
//{
//	std::uniform_int_distribution<int> distribution(0,1);

//  for(int i = 0; i < constants::V; i++)
//  {
//    sig_link[i][0] = 2*distribution(generator1)-1;// -1 or 1
//    sig_link[i][1] = 2*distribution(generator1)-1;
//  }
//	return;
//}



