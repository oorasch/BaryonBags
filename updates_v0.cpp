#include<iostream>
#include<algorithm>
#include<vector>
#include<gsl/gsl_linalg.h>
#include "parameters.h"
#include "updates.h"

using namespace std;

//routine for one sweep
void do_updates(const int neib[][4], )
{
	for(int x = 0; x < param::V; x++)
	{
		int y = neib[x][1];
		vector<int> x_links = {dimer.at(x).at(0), dimer.at(x).at(1), dimer.at(neib[x][2]).at(0), dimer.at(neib[x][3]).at(1)},
								y_links = {dimer.at(y).at(0), dimer.at(y).at(1), dimer.at(neib[y][2]).at(0), dimer.at(neib[y][3]).at(1)};
		
		update_routine(x, x_prime, monomer.at(x), monomer.at(y), x_links, y_links);
		
		y = neib[x][0];
		update_routine(x, x_prime, monomer.at(x), monomer.at(y), x_links, y_links);
		
	}
}

void update_routine(const int& x, const int& y, int& mx, int& my, int& d, vector<int>& x_links, vector<int>& y_links)
{
	int site_type_x = check_site_type(mx, x_links),
			site_type_y = check_site_type(my, y_links);
		
	if(site_type_x == 0 && site_type_y == 0 && d == 1)
	{
		//change to: 4, 4, 0
	}
	else if(site_type_x == 4 && site_type_y == 4 && d == 0)
	{
		//change to: 0, 0, 1
	}
	else if(site_type_x == 0 && site_type_y == 0 && d == 2 )
	{
		//change to: 1, 1, 1
	}
	else if(site_type_x == 1 && site_type_y == 1 && d == 1)
	{
		//change to: 0, 0, 2
	}
	else if(site_type_x == 1 && site_type_y == 0 && d == 1)
	{
		//change to: 2, 4, 0
	}
	else if(site_type_x == 2 && site_type_y == 4 && d == 0)
	{
		//change to: 1, 0, 1
	}
	else if(site_type_x == 2 && site_type_y == 0 && d == 1)
	{
		//change to: 3 (det!), 4, 0
	}
	else if(site_type_x == 3 && site_type_y == 4 && d == 0)
	{
		//change to: 2 (det!), 0, 1
	}
	else if(site_type_x == 4 && site_type_y == 0 && d == 2)
	{
		//change to: 2, 1, 1
	}
	else if(site_type_x == 2 && site_type_y == 1 && d == 1)
	{
		//change to: 4, 0, 2
	}
	else if(site_type_x == 5 && site_type_y == 0 && d == 1)
	{
		//change to: 1, 4, 0
	}
	else if(site_type_x == 1 && site_type_y == 4 && d == 0)
	{
		//change to: 5, 0, 1
	}
	else if(site_type_x == 1 && site_type_y == 1 && d == 0)
	{
		//change to: 5, 5, 1
	}
	else if(site_type_x == 5 && site_type_y == 5 && d == 1)
	{
		//change to: 1, 1, 0
	}
	else if(site_type_x == 2 && site_type_y == 1 && d == 0)
	{
		//change to: 1, 5, 1
	}
	else if(site_type_x == 1 && site_type_y == 5 && d == 1)
	{
		//change to: 2, 1, 0
	}
	else if(site_type_x == 3 && site_type_y == 1 && d == 0)
	{
		//change to: 2 (det!), 5, 1
	}
	else if(site_type_x == 2 && site_type_y == 5 && d == 1)
	{
		//change to: 3 (det!), 1, 0
	}
	else if(site_type_x == 2 && site_type_y == 2 && d == 1)
	{
		//change to: 3, 3, 0 (det!)
	}
	else if(site_type_x == 3 && site_type_y == 3 && d == 0)
	{
		//change to: 2, 2, 1 (det!)
	}
	else if(site_type_x == 4 && site_type_y == 4 && d == 2)
	{
		//change to: 6, 6, 3 (det!)
	}
	else if(site_type_x == 6 && site_type_y == 6 && d == 3)
	{
		//change to: 4, 4, 2 (det!)
	}
	else if(site_type_x == 4 && site_type_y == 4 && d == 0)
	{
		//change to: 0, 0, 1
	}
	else if(site_type_x == 0 && site_type_y == 0 && d == 1)
	{
		//change to: 4, 4, 0
	}
	else if(site_type_x == 1 && site_type_y == 1 && d == 1)
	{
		//change to: 0, 0, 2
	}
	else if(site_type_x == 0 && site_type_y == 0 && d == 2)
	{
		//change to: 1, 1, 1
	}
	else if(site_type_x == 4 && site_type_y == 2 && d == 0)
	{
		//change to: 0, 0, 1
	}
	else if(site_type_x == 4 && site_type_y == 3 && d == 0)
	{
		//change to: 0, 2 (det!), 1
	}
	else if(site_type_x == 0 && site_type_y == 2 && d == 1)
	{
		//change to: 4, 3, 0 (det!)
	}
	else if(site_type_x == 4 && site_type_y == 1 && d == 0)
	{
		//change to: 0, 5, 1
	}
	else if(site_type_x == 0 && site_type_y == 5 && d == 1)
	{
		//change to: 4, 1, 0
	}
	else if(site_type_x == 5 && site_type_y == 5 && d == 1)
	{
		//change to: 1, 1, 0
	}
	else if(site_type_x == 1 && site_type_y == 1 && d == 0)
	{
		//change to: 
	}
	else if(site_type_x == 5 && site_type_y == 1 && d == 1)
	{
		//change to: 1, 2, 0
	}
	else if(site_type_x == 1 && site_type_y == 2 && d == 0)
	{
		//change to: 5, 1, 1
	}
	else if(site_type_x == 5 && site_type_y == 2 && d == 1)
	{
		//change to: 1, 3, 0
	}
	else if(site_type_x == 1 && site_type_y == 3 && d == 0)
	{
		//change to: 5, 2, 1
	}
	else if(site_type_x == 3 && site_type_y == 3 && d == 0)
	{
		//change to: 2, 2, 1 (det!)
	}
	else if(site_type_x == 2 && site_type_y == 2 && d == 1)
	{
		//change to: 3, 3, 0 (det!)
	}
	else if(site_type_x == 6 && site_type_y == 6 && d == 3)
	{
		//change to: 4, 4, 2 (det!)
	}
	else if(site_type_x == 4 && site_type_y == 4 && d == 2)
	{
		//change to: 6, 6, 3 (det!)
	}
	else if(site_type_x == 3 && site_type_y == 2 && d == 0)
	{
		//change to: 2 (det!), 1, 1
	}
	else if(site_type_x == 2 && site_type_y == 1 && d == 1)
	{
		//change to: 3 (det!), 2, 0
	}
	else if(site_type_x == 2 && site_type_y == 1 && d == 1)
	{
		//change to: 3 (det!), 2, 0
	}
	else if(site_type_x == && site_type_y == && d == )
	{
		//change to: 
	}


}

int check_site_types(const int& m, const vector<int>& links)
{
	int site type = -1;
	
	if(m == 0)
	{
		if(links.at(0) == 3 || links.at(1) == 3 || links.at(2) == 3 || links.at(3) == 3) site_type = 6;
		else if(links.at(0) == 2 || links.at(1) == 2 || links.at(2) == 2 || links.at(3) == 2) site_type = 0;
		else site_type = 5;
	}
	else if(m == 1)
	{
		if(links.at(0) == 2 || links.at(1) == 2 || links.at(2) == 2 || links.at(3) == 2) site_type = 4;
		else site_type = 1;
	}
	else if(m == 2)
	{
		site_type = 2;
	}
	else site_type = 3;
	
	if(site_type < 0)
	{
		abort();
	}
	
	return site_type;
}
