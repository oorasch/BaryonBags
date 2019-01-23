#include<iostream>
#include<algorithm>
#include<vector>
#include<gsl/gsl_linalg.h>
#include "parameters.h"
#include "measure.h"

using namespace std;

double meas_chir_cond(const vector<int>& site_type)
{
	vector<int> monomers = {0,1,2,1,0,0};
	double tmp = 0.0;
	
	for(int i = 0; i < param::V; i++)
	{
		tmp += monomers.at(site_type.at(i))/(2.0*param::m);
//		if(site_type.at(i) == 5) tmp += 3.0/2.0/param::m;
	}
	
	return tmp;
}

double meas_chir_susc(const vector<int>& site_type)
{
	vector<int> monomers = {0,1,2,1,0,0};
	double tmp = 0.0;
	
	for(int i = 0; i < param::V; i++)
	{
		tmp += monomers.at(site_type.at(i));
	}
	
	return (tmp*tmp/param::m/param::m) - tmp*(tmp-1.0)/param::m/param::m;
	
}
