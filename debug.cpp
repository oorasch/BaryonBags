#include<iostream>
#include<algorithm>
#include<vector>
#include<gsl/gsl_linalg.h>
#include "parameters.h"
#include "debug.h"

using namespace std;

void print_md_occ(const vector<int>& site_type, const vector<vector<int>>& dimer)
{
	for(int i = 0; i < param::V; i++) cout << "@ x = " << i << ": (" << site_type.at(i) << "," << dimer.at(i).at(0) << "," << dimer.at(i).at(1) << ")" << endl;
}

void print_lattice_occ(const vector<int>& site_type, const vector<vector<int>>& dimer)
{
	for(int i = 0; i < param::V; i++)
	{
		 cout <<  site_type.at(i) << " ";
		 if((i+1)%param::Ns == 0) cout << endl;
	}
}

void check_occupation (const int neib[][4], const vector<int>& site_type, const vector<vector<int>>& dimer)
{
	vector<int> monomers = {0,1,2,1,0,0};
	
	for(int i = 0; i < param::V; i++)
	{
		int tmp = monomers.at(site_type.at(i));
		tmp += dimer.at(i).at(0) + dimer.at(i).at(1);
		tmp += dimer.at(neib[i][2]).at(0) + dimer.at(neib[i][3]).at(1);

		
		if(tmp%3 != 0)
		{
			cout << "Fail @ x = " << i << ": " << monomers.at(site_type.at(i)) << ", " << dimer.at(i).at(0);
			cout << ", " << dimer.at(i).at(1) << ", " << dimer.at(neib[i][2]).at(0) << ", " << dimer.at(neib[i][3]).at(1) << endl;
		}
	}

}


void fill_bag_vec(const int neib[][4], const vector<int>& monomer, const vector<vector<int>>& dimer, vector<vector<int>>& bag_vec)
{
	bag_vec.clear();
	
	for(int x = 0; x < param::V; x++)
	{
		if(monomer.at(x) == 0 && dimer.at(x).at(0) == 0 && dimer.at(x).at(1) == 0)
		{
			if(bag_vec.empty()) bag_vec.push_back(vector<int>(1, x));
			else{
			
				bool new_bag = true;
								
				for(size_t k = 0; k < bag_vec.size(); k++)
				{
					if(is_neib(neib, bag_vec.at(k), x))
					{
						bag_vec.at(k).push_back(x);
						sort(bag_vec.at(k).begin(), bag_vec.at(k).end());
						new_bag = false;
						break;
					}
				}
				
				if(new_bag) bag_vec.push_back(vector<int>(1, x));
				
			}
		}
	}

	return;
}

//check if x is a neighbor to any site in the bag
bool is_neib(const int neib[][4], const vector<int>& bag_coords, const int& x)
{
	bool is_neib = false;

	for(auto it1: bag_coords)
	{	
		int neib0 =  neib[it1][0], neib1 =  neib[it1][1], neib2 =  neib[it1][2], neib3 =  neib[it1][3];
		
		if(x == neib0)
		{
			is_neib = true;
			break;
		}
		else if(x == neib1)
		{
			is_neib = true;
			break;
		}
		else if(x == neib2)
		{
			is_neib = true;
			break;
		}
		else if(x == neib3)
		{
			is_neib = true;
			break;
		}
			
//		if(!is_neib) cout << "Alert: disconnected bag! (" << it1 << ")" << endl;
	}
	
	return is_neib;
}



