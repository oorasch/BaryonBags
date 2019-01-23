#include<iostream>
#include<algorithm>
#include<random>
#include<chrono>
#include<vector>
#include<gsl/gsl_linalg.h>
#include "parameters.h"
#include "updates.h"
#include "debug.h"
#include "bags.h"

using namespace std;

const unsigned seed(std::chrono::system_clock::now().time_since_epoch().count());
std::default_random_engine generator(seed);

//routine for one sweep, change locally site_type and dimer occupation
void do_updates(const int& sweeps, const int neib[][4], vector<int>& site_type, vector<vector<int>>& dimer, vector<int>& bag_vec, vector<vector<bool>>& bag_sites, vector<double>& bag_det)
{
	std::uniform_int_distribution<int> dist1(0, param::V-1);
	std::uniform_int_distribution<int> dist2(0,1);
	
	int x = 0, y = 0, dir = 0;

	for(int n = 0; n < sweeps; n++)
	{
		for(int i = 0; i < param::V; i++)
		{			
			x = dist1(generator);
			dir = dist2(generator);
		  y = neib[x][dir];
		  
			//site_type.at(x), site_type.at(y), dimer.at(x).at(1)
			B_bar_update(site_type.at(x), site_type.at(y), dimer.at(x).at(dir));
			
			x = dist1(generator), dir = dist2(generator);
			y = neib[x][dir];
			B_update(x, y, neib, site_type.at(x), site_type.at(y), dimer.at(x).at(dir), bag_vec, bag_sites, bag_det);
			
			if(param::max_label != 0) cout << param::max_label << endl;
			
//			print_lattice_occ(site_type, dimer);
			
		}
		
	}
}

//here we update only the complimentary region
void B_bar_update(int& site_type_x, int& site_type_y, int& d)
{
	double rho = 0.0;
	double mm2 = 2.0*param::m;
	std::uniform_real_distribution<double> distribution(0,1);
	double rr = distribution(generator);//random number
	
//	cout << rr << endl;
	
	if(site_type_x == 0 && site_type_y == 0 && d == 1)
	{
		//change to: 3, 3, 0
		rho = 3.0*mm2*mm2;
		
		if(rr < rho)
		{
			site_type_x = 3;
			site_type_y = 3;
			d = 0;
		}
		
	}
	else if(site_type_x == 3 && site_type_y == 3 && d == 0)
	{
		//change to: 0, 0, 1
		rho = 1.0/(3.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 0;
			site_type_y = 0;
			d = 1;
		}
	}
	else if(site_type_x == 0 && site_type_y == 0 && d == 2 )
	{
		//change to: 1, 1, 1
		rho = 4.0*mm2*mm2;
		
		if(rr < rho)
		{
			site_type_x = 1;
			site_type_y = 1;
			d = 1;
		}
	}
	else if(site_type_x == 0 && site_type_y == 1 && d == 1)
	{
		//change to: 3, 2, 0
		rho = 3.0*mm2*mm2/2.0;
		
		if(rr < rho)
		{
			site_type_x = 3;
			site_type_y = 2;
			d = 0;
		}
	}
	else if(site_type_x == 3 && site_type_y == 2 && d == 0)
	{
		//change to: 0, 1, 1
		rho = 2.0/(3.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 0;
			site_type_y = 1;
			d = 1;
		}
	}
	else if(site_type_x == 2 && site_type_y == 3 && d == 0)
	{
		//change to: 1, 0, 1
		rho = 2.0/(3.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 1;
			site_type_y = 0;
			d = 1;
		}
	}
	else if(site_type_x == 1 && site_type_y == 0 && d == 1)
	{
		//change to: 2, 3, 0
		rho =  3.0*mm2*mm2/2.0;
		
		if(rr < rho)
		{
			site_type_x = 2;
			site_type_y = 3;
			d = 0;
		}
	}
	else if(site_type_x == 0 && site_type_y == 4 && d == 1)
	{
		//change to: 3, 1, 0
		rho = 3.0*mm2*mm2;
		
		if(rr < rho)
		{
			site_type_x = 3;
			site_type_y = 1;
			d = 0;
		}
	}
	else if(site_type_x == 3 && site_type_y == 1 && d == 0)
	{
		//change to: 0, 4, 1
		rho = 1.0/(3.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 0;
			site_type_y = 4;
			d = 1;
		}
	}
	else if(site_type_x == 1 && site_type_y == 3 && d == 0)
	{
		//change to: 4, 0, 1
		rho = 1.0/(3.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 4;
			site_type_y = 0;
			d = 1;
		}
	}
	else if(site_type_x == 4 && site_type_y == 0 && d == 1)
	{
		//change to: 1, 3, 0
		rho = 3.0*mm2*mm2;
		
		if(rr < rho)
		{
			site_type_x = 1;
			site_type_y = 3;
			d = 0;
		}
	}
	else if(site_type_x == 1 && site_type_y == 1 && d == 0)
	{
		//change to: 4, 4, 1
		rho = 1.0/(3.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 4;
			site_type_y = 4;
			d = 1;
		}
	}
	else if(site_type_x == 4 && site_type_y == 4 && d == 1)
	{
		//change to: 1, 1, 0
		rho = 3.0*mm2*mm2;
		
		if(rr < rho)
		{
			site_type_x = 1;
			site_type_y = 1;
			d = 0;
		}
	}
	else if(site_type_x == 1 && site_type_y == 1 && d == 1)
	{
		//change to: 0, 0, 2 or 2, 2, 0 with p = 0.5
		std::uniform_real_distribution<double> distribution(0,1);
		double rr_tmp = distribution(generator);
		
		if(rr_tmp < 0.5)
		{
			rho = 1.0/(4.0*mm2*mm2);
			
			if(rr < rho)
			{
				site_type_x = 0;
				site_type_y = 0;
				d = 2;
			}
		}
		else{
			
			rho = 3.0*mm2*mm2/4.0;
		
			if(rr < rho)
			{
				site_type_x = 2;
				site_type_y = 2;
				d = 0;
			}
			
		}
	}
	else if(site_type_x == 2 && site_type_y == 2 && d == 0)
	{
		//change to: 1, 1, 1
		rho = 4.0/(3.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 1;
			site_type_y = 1;
			d = 1;
		}
		
	}
	else if(site_type_x == 1 && site_type_y == 2 && d == 0)
	{
		//change to: 4, 1, 1
		rho = 2.0/(3.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 4;
			site_type_y = 1;
			d = 1;
		}
	}
	else if(site_type_x == 4 && site_type_y == 1 && d == 1)
	{
		//change to: 1, 2, 0
		rho = 3.0*mm2*mm2/2.0;
		
		if(rr < rho)
		{
			site_type_x = 1;
			site_type_y = 2;
			d = 0;
		}
	}
	else if(site_type_x == 1 && site_type_y == 4 && d == 1)
	{
		//change to: 2, 1, 0
		rho = 3.0*mm2*mm2/2.0;
		
		if(rr < rho)
		{
			site_type_x = 2;
			site_type_y = 1;
			d = 0;
		}
	}
	else if(site_type_x == 2 && site_type_y == 1 && d == 0)
	{
		//change to: 1, 4, 1
		rho = 2.0/(3.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 1;
			site_type_y = 4;
			d = 1;
		}
	}
	else if(site_type_x == 0 && site_type_y == 3 && d == 2)
	{
		//change to: 1, 2, 1
		rho = 2.0*mm2*mm2;
		
		if(rr < rho)
		{
			site_type_x = 1;
			site_type_y = 2;
			d = 1;
		}
	}
	else if(site_type_x == 1 && site_type_y == 2 && d == 1)
	{
		//change to: 0, 3, 2
		rho = 1.0/(2.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 0;
			site_type_y = 3;
			d = 2;
		}
	}
	else if(site_type_x == 2 && site_type_y == 1 && d == 1)
	{
		//change to: 3, 0, 2
		rho = 1.0/(2.0*mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 3;
			site_type_y = 0;
			d = 2;
		}
	}
	else if(site_type_x == 3 && site_type_y == 0 && d == 2)
	{
		//change to: 2, 1, 1
		rho = 2.0*mm2*mm2;
		
		if(rr < rho)
		{
			site_type_x = 2;
			site_type_y = 1;
			d = 1;
		}
	}
	else if(site_type_x == 2 && site_type_y == 2 && d == 1)
	{
		//change to: 3, 3, 2
		rho = 1.0/(mm2*mm2);
		
		if(rr < rho)
		{
			site_type_x = 3;
			site_type_y = 3;
			d = 2;
		}
	}
	else if(site_type_x == 3 && site_type_y == 3 && d == 2)
	{
		//change to: 2, 2, 1
		rho = mm2*mm2;
		
		if(rr < rho)
		{
			site_type_x = 2;
			site_type_y = 2;
			d = 1;
		}
	}
	
}

//here we update only the bag region
void B_update(const int& x, const int& y, const int neib[][4], int& site_type_x, int& site_type_y, int& d, vector<int>& bag_vec, vector<vector<bool>>& bag_sites, vector<double>& bag_det)
{
	int site_type_x_prime = 0, site_type_y_prime = 0, d_prime = 0;
	double mm2 = 2.0*param::m;
	double rho = 0.0;
	std::uniform_real_distribution<double> distribution(0,1);
	double rr = distribution(generator);//random number
	vector<vector<int>> local_bags; //temp storage of the bag neighborhood of a site, max length = 4, first index neibs, 0 = site, 1 = bag number 
//	vector<int> 
	
	if(site_type_x == 0 && site_type_y == 2 && d == 1)
	{
		//change to: 3, 5, 0
		site_type_x_prime = 3;
		site_type_y_prime = 5;
		d_prime = 0;
		rho = 1.0/mm2; // comes from change in compl. region
		//add site y to bag, check neigbors
		//check neighbors
		
		check_local_bag_occupation(y, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //make a new bag
		{
			rho *= mm2*mm2*mm2;
			
			cout << param::m  << " " << rho << endl;
			
			if(rr < rho)
			{
				make_bag(y, bag_vec, bag_sites, bag_det);
				
				site_type_x = 3;
				site_type_y = 5;
				d = 0;
			}
			
		}
//		if else (local_bags.size() == 1)//add site
//		{
//			//bag number = local_bags.at(0).at(1)
//			//neib site = local_bags.at(0).at(0)
//			bag_sites_prime = bag_sites.at(local_bags.at(0).at(1))
//			bag_site_prime.at(y) = 1;
//			
//			rho *=;// dets
//			
//			if(rr < rho)
//			{
//				bag_size(local_bags.at(0).at(1))++;
//				site_type_x = 3;
//				site_type_y = 5;
//				d = 0;
//			}
//			
//		}
//		if else (local_bags.size() == 2)
//		{
//		
//		}
//		if else (local_bags.size() == 3)
//		{
//		
//		}


//		
//		
	}
	else if(site_type_x == 3 && site_type_y == 5 && d == 0)
	{
		//change to: 0, 2, 1
		site_type_x_prime = 0;
		site_type_y_prime = 2;
		d_prime = 1;
		//remove site y from bag
		rho = mm2;// comes from change in compl. region
		
		check_local_bag_occupation(y, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //delete bag
		{
			rho /= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				del_bag_1(y, bag_vec, bag_sites, bag_det);
				
				site_type_x = 0;
				site_type_y = 2;
				d = 1;
			}
			
		}
//		if else (local_bags.size() == 1)//add site
//		{
//			//bag number = local_bags.at(0).at(1)
//			//neib site = local_bags.at(0).at(0)
//			bag_sites_prime = bag_sites.at(local_bags.at(0).at(1))
//			bag_site_prime.at(y) = 1;
//			
//			rho *=;// dets
//			
//			if(rr < rho)
//			{
//				bag_size(local_bags.at(0).at(1))++;
//				site_type_x = 3;
//				site_type_y = 5;
//				d = 0;
//			}
//			
//		}
//		if else (local_bags.size() == 2)
//		{
//		
//		}
//		if else (local_bags.size() == 3)
//		{
//		
//		}
	}
	else if(site_type_x == 2 && site_type_y == 0 && d == 1)
	{
		//change to: 5, 3, 0
		site_type_x_prime = 5;
		site_type_y_prime = 3;
		d_prime = 0;
		//add site x to bag
		rho = 1.0/mm2;// comes from change in compl. region
		
//		cout << "Hi x = " << x << " " << param::max_label << endl;
		check_local_bag_occupation(x, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //make new bag
		{
			rho *= mm2*mm2*mm2;
			
			if(rr < rho)
			{
//			  cout << "Hi again x = " << x << " " << param::max_label <<endl;
				make_bag(x, bag_vec, bag_sites, bag_det);
//				 cout << "Hi again 2 x = " << x << " " << param::max_label << endl;
				
				site_type_x = 5;
				site_type_y = 3;
				d = 0;
			}
			
		}
	}
	else if(site_type_x == 5 && site_type_y == 3 && d == 0)
	{
		//change to: 2, 0 1
		site_type_x_prime = 2;
		site_type_y_prime = 0;
		d_prime = 1;
		//remove site x from bag
		rho = mm2;// comes from change in compl. region
		
		check_local_bag_occupation(x, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //delete bag
		{
			rho /= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				del_bag_1(x, bag_vec, bag_sites, bag_det);
				
				site_type_x = 2;
				site_type_y = 0;
				d = 1;
			}
			
		}
	}
	else if(site_type_x == 1 && site_type_y == 2 && d == 1)
	{
		//change to: 2, 5, 0
		site_type_x_prime = 2;
		site_type_y_prime = 5;
		d_prime = 0;
		//add site y to bag
		rho = 1.0/(2.0*mm2);// comes from change in compl. region
		
		check_local_bag_occupation(y, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //make new bag
		{
			rho *= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				make_bag(y, bag_vec, bag_sites, bag_det);
				
				site_type_x = 2;
				site_type_y = 5;
				d = 0;
			}
			
		}
	}
	else if(site_type_x == 2 && site_type_y == 5 && d == 0)
	{
		//change to: 1, 2, 1
		site_type_x_prime = 1;
		site_type_y_prime = 2;
		d_prime = 1;
		//remove site y from bag
		rho = 2.0*mm2;// comes from change in compl. region
		
		check_local_bag_occupation(y, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //delete bag
		{
			rho /= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				del_bag_1(y, bag_vec, bag_sites, bag_det);
				
				site_type_x = 1;
				site_type_y = 2;
				d = 1;
			}
			
		}
		
	}
	else if(site_type_x == 2 && site_type_y == 1 && d == 1)
	{
		//change to: 5, 2, 0
		site_type_x_prime = 5;
		site_type_y_prime = 2;
		d_prime = 0;
		//add site x to bag
		rho = 1.0/(2.0*mm2);// comes from change in compl. region
		
		check_local_bag_occupation(x, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //make new bag
		{
			rho *= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				make_bag(x, bag_vec, bag_sites, bag_det);
				
				site_type_x = 5;
				site_type_y = 2;
				d = 0;
			}
			
		}
	}
	else if(site_type_x == 5 && site_type_y == 2 && d == 0)
	{
		//change to: 2, 1, 1
		site_type_x_prime = 2;
		site_type_y_prime = 1;
		d_prime = 1;
		//remove site x from bag
		rho = 2.0*mm2;// comes from change in compl. region
		
		check_local_bag_occupation(x, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //delete bag
		{
			rho /= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				del_bag_1(x, bag_vec, bag_sites, bag_det);
				
				site_type_x = 2;
				site_type_y = 1;
				d = 1;
			}
			
		}
		
	}
	else if(site_type_x == 2 && site_type_y == 2 && d == 1)
	{
//		//change to: 5, 5, 0
//		site_type_x_prime = 5;
//		site_type_y_prime = 5;
//		d_prime = 0;
//		//add site x and y to bag
//		rho = 1.0/9.0/(4.0*param::m*param::m);// comes from change in compl. region
	}
	else if(site_type_x == 5 && site_type_y == 5 && d == 0)
	{
//		//change to: 2, 2, 1
//		site_type_x_prime = 2;
//		site_type_y_prime = 2;
//		d_prime = 1;
//		//remove site x and y from bag
//		rho = 9.0*(4.0*param::m*param::m);// comes from change in compl. region
	}
	else if(site_type_x == 4 && site_type_y == 2 && d == 1)
	{
		//change to: 1, 5, 0
		site_type_x_prime = 1;
		site_type_y_prime = 5;
		d_prime = 0;
		//add site y to bag
		rho = 1.0/mm2;// comes from change in compl. region

		check_local_bag_occupation(y, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //
		{
			rho *= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				make_bag(y, bag_vec, bag_sites, bag_det);
				
				site_type_x = 1;
				site_type_y = 5;
				d = 0;
			}
		}
		
	}
	else if(site_type_x == 1 && site_type_y == 5 && d == 0)
	{
		//change to: 4, 2, 1
		site_type_x_prime = 4;
		site_type_y_prime = 2;
		d_prime = 1;
		//add site y to bag
		rho = mm2;// comes from change in compl. region
		
		check_local_bag_occupation(y, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //
		{
			rho /= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				del_bag_1(y, bag_vec, bag_sites, bag_det);
				
				site_type_x = 4;
				site_type_y = 2;
				d = 1;
			}
		}
		
	}
	else if(site_type_x == 2 && site_type_y == 4 && d == 1)
	{
		//change to: 5, 1, 0
		site_type_x_prime = 5;
		site_type_y_prime = 1;
		d_prime = 0;
		//add site x to bag
		rho = 1.0/mm2;// comes from change in compl. region

		check_local_bag_occupation(x, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //
		{
			rho *= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				make_bag(x, bag_vec, bag_sites, bag_det);
				
				site_type_x = 5;
				site_type_y = 1;
				d = 0;
			}
		}
	}
	else if(site_type_x == 5 && site_type_y == 1 && d == 0)
	{
		//change to: 2, 4, 1
		site_type_x_prime = 2;
		site_type_y_prime = 4;
		d_prime = 1;
		//remove site x from bag
		rho = mm2;// comes from change in compl. region
		
		check_local_bag_occupation(x, neib, bag_vec, local_bags);
		
		if(local_bags.size() == 0) //
		{
			rho /= mm2*mm2*mm2;
			
			if(rr < rho)
			{
				del_bag_1(x, bag_vec, bag_sites, bag_det);
				
				site_type_x = 2;
				site_type_y = 4;
				d = 1;
			}
		}
	}
	
}

//check how many neighbourign sites are occupied by bags
void check_local_bag_occupation(const int& site, const int neib[][4], const vector<int>& bag_vec, vector<vector<int>>& local_bags)
{
	vector<int> tmp (2, 0);
	local_bags.clear();
	for(int i = 0; i < 4; i++)
	{
		if(bag_vec.at(neib[site][i]) != 0)
		{
			tmp.at(0) = site; //store neighboring site
			tmp.at(1) = bag_vec.at(neib[site][i]);//store bag num
			
			local_bags.push_back(tmp);
		}
	}
}

////check how many neighbors of site n belong to a bag (not necessarily the same)
//bool check_connection(const int& nn, const int& mm, const vector<vector<bool>>& bag_site_vec)
//{

//	int qq = -1;
//	
//	if(bag_site_vec.at(nn).at(0)) qq = neib[nn][0];
//	 	else qq = nn;
//	 		
//	while(qq != nn || qq != mm)
//	{
//		for(int j = 0; j < 4; j++)
//		{
//	 		if(bag_site_vec.at(nn).at(j)) qq = neib[nn][j];
//		}
//		
//	}
//	

//}


