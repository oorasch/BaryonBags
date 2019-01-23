#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<random>
//#include "bags.h"
#include "updates.h"
#include "parameters.h"
#include "initialize.h"
#include "measure.h"
#include "debug.h"
#include "bags.h"

using namespace std;

int main(){

vector<int> monomer, site_type (param::V, -1);
vector<vector<int>>  dimer (param::V, vector<int>(2, 0));
vector<int> bag_vec (param::V, 0); //stores the labels of bags, 0 = not in a bag, l_bag = 1, 2, ..., V, max number of bags = V
vector<int> bag_size (param::V, 0); //stores the size of bags, ib = 
vector<vector<bool>> bag_sites (param::V, vector<bool>(param::V, false)); //first index is bag number, second sites in bag, 0 means not in this bag, or in no bag at all
vector<double> bag_det (param::V, 0); //stores the determinant of bags
int neib [param::V][4];

neib_init(param::Nt, param::Ns, neib);
init_md_config(neib, site_type, dimer, bag_vec);

//print_lattice_occ(site_type, dimer);

//ofstream outfile1, outfile2;
//outfile1.open("cond_data.dat");
//outfile2.open("susc_data.dat");

//for(int i = 0; i < 10; i++)
//{
//	double tmp1 = 0.0, tmp2 = 0.0;
//	param::m += 1.0/(10.0);

//	do_updates(param::nskip, neib, site_type, dimer, bag_vec, bag_sites, bag_det);

//	for(int i = 0; i < param::nmeas; i++)
//	{
//		do_updates(param::nskip, neib, site_type, dimer, bag_vec, bag_sites, bag_det);
//		tmp1 += meas_chir_cond(site_type);
//		tmp2 += meas_chir_susc(site_type);
//	}
//	
//	tmp1 /= 1.0*param::nmeas;
//	tmp2 /= 1.0*param::nmeas;
//	
//	outfile1 << param::m << " " << tmp1 << endl;
//	outfile2 << param::m << " " << tmp2 << endl;

//	check_occupation(neib, site_type, dimer);

////	print_md_occ(site_type, dimer);
//}
//outfile1.close(); outfile2.close();

vector<bool> test (param::V, false);

test.at(1) = true;
test.at(5) = true;
test.at(13) = true;
//test.at(8) = true;


cout << calc_det(neib, test) << endl;


return 0;
}
