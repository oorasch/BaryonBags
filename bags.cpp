#include<iostream>
#include<algorithm>
#include<vector>
#include<gsl/gsl_linalg.h>
#include "parameters.h"
#include "bags.h"


using namespace std;


void make_bag(const int& site, vector<int>& bag_vec, vector<vector<bool>>& bag_sites, vector<double>& bag_det)
{
	param::max_label++;
	bag_vec.at(site) = param::max_label;
	bag_sites.at(param::max_label).at(site) = true;
	bag_det.at(param::max_label) = 8.0*param::m*param::m*param::m;
}

//deletes a 1-site bag
void del_bag_1(const int& site, vector<int>& bag_vec, vector<vector<bool>>& bag_sites, vector<double>& bag_det)
{
	int label = bag_vec.at(site); 
	param::max_label--;
	bag_vec.at(site) = 0;
	bag_sites.at(label).at(site) = false;
	bag_det.at(label) = 1.0;
	
	//We delete the the bag with label from bag_sites and add a new empty element at the back of the vector
	//This should only be a temporary solution. Better would be to use a label_vec that stores the canonical labels for the bags
	bag_sites.erase(bag_sites.begin()+label);
	bag_sites.push_back(vector<bool>(param::V, false));
	
	//reorder all bags: 1,2,3,5,6, ... -> 1,2,3,4,5,...
	for(int i = 0; i < param::V; i++)
	{
		if(bag_vec.at(i) > i) bag_vec.at(i)--;
	}
	
}


void add_1_site(const int& site, const int& bag_label, vector<bool>& bag_sites_bag, vector<double>& bag_det)
{
	bag_sites_bag.at(site) = true;
	bag_det.at(bag_label) = 0.0; // recalculate the bag determinant

}

void rem_1_site(const int& site, const int& bag_label, vector<bool>& bag_sites_bag, vector<double>& bag_det)
{
	bag_sites_bag.at(site) = false;
	bag_det.at(bag_label) = 0.0; // recalculate the bag determinant

}

double calc_det(const int neib [][4], const vector<bool>& bag)
{
	int k = 0, gamma = 0;
	double bag_det, mm2 = 2.0*param::m;
	vector<double> D;
	
	for(int kx = 0; kx < param::V; kx++)
	{
		for(int ky = 0; ky < param::V; ky++)
		{
//			int k = kx + ky*param::V;
			
			if(bag.at(kx) && bag.at(ky))
			{
				if((kx%param::Ns)%2 == 0) gamma = 1;
				else gamma = -1;
				
				if(kx == ky) D.push_back(mm2*mm2*mm2);
				else if(ky == neib[kx][0]) D.push_back(gamma*1.0); //+ from kernel
				else if(ky == neib[kx][1]) D.push_back(1.0); //+ from kernel, + from gamma
				else if(kx == neib[ky][0]) D.push_back(gamma*(-1.0)); //- from kernel
				else if(kx == neib[ky][1]) D.push_back(-1.0); //- from kernel, + from gamma
				else D.push_back(0.0);
			}
		}
	}
	
	int tmp=0, bag_size = sqrt(D.size());
	cout << bag_size << endl;
	for(auto it: D)
	{
	 	cout << it << " ";
	 	if((tmp+1)%bag_size == 0) cout << endl;
	 	tmp++;
	}
	
	if(bag_size == 1) bag_det = mm2*mm2*mm2;
	else if(bag_size == 2) bag_det = (mm2*mm2*mm2)*(mm2*mm2*mm2) + 1.0;
	else if (bag_size == 3) bag_det = (mm2*mm2*mm2)*(mm2*mm2*mm2)*(mm2*mm2*mm2) + 2.0*(mm2*mm2*mm2)*(mm2*mm2*mm2);
	else
	{
		int s;
		gsl_matrix_view m = gsl_matrix_view_array (D.data(), bag_size, bag_size);
		gsl_permutation * p = gsl_permutation_alloc (sqrt(D.size()));
		gsl_linalg_LU_decomp (&m.matrix, p, &s);
		
		bag_det = gsl_linalg_LU_det(&m.matrix, s);
		gsl_permutation_free(p);
	}
  
  return bag_det;
 
	
}

//void add_2_site(const int& site, const int& bag_label_1, const int& bag_label_2, vector<vector<bool>>& bag_sites, vector<double>& bag_det)
//{
//	int fin_bag_label = 0, // final bag label
//		  del_bag_label = 0; // bag that gets absorbed
//	
//	if (bag_label_1 == bag_label_2)//sites belong to same bag
//	{
//		fin_bag_label = bag_label_1;
//	}
//	else if (bag_label_1 < bag_label_2)
//	{
//		fin_bag_label = bag_label_1;
//		del_bag_label = bag_label_2;
//	}
//	else //if (bag_label_2 > bag_label_1)
//	{
//		fin_bag_label = bag_label_2;
//		del_bag_label = bag_label_1;
//	}
//	
//	if(del_bag_label == 0)//sites gets added to a bulky bag
//	{
//		bag_sites.at(fin_bag_label).at(site) = true;
//		bag_det.at(bag_label) = 0.0; // recalculate the bag determinant
//	}
//	else // site merges two previously separated bags
//	{
//		//add site to bag
//		bag_sites.at(fin_bag_label).at(site) = true;
//		
//		//copy all bag sites to the one with smaller label
//		for(int i = 0; i < param::V; i++)
//		{
//			bag_sites.at(fin_bag_label).at(i) = bag_sites.at(del_bag_label).at(i);
//		}
//		
//		bag_det.at(bag_label) = 0.0; // recalculate the bag determinant
//		
//		//delete bag
//		bag_sites.erase(bag_sites.begin()+label);
//		bag_sites.push_back(vector<bool>(param::V, false));
//	}

//}







