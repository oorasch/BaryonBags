#include<iostream>
#include<algorithm>
#include<vector>
#include<gsl/gsl_linalg.h>

using namespace std;

const int Nt = 4, Ns = 4;
const int V = Nt*Ns;
const double m = 1.0;


class bag
{
	vector<int> bag_coords;
	double det_D;
	
	public:
		bag(); //default constructor
		~bag(); // default destructor
		void set_bag_coords(vector<int> vec); // set the bag coordinates
		vector<int> get_bag_coords(); // get the bag coordinates
		void set_det_D(double det_M); // set the bag determinant
		double get_det_D(); // get the bag determinant
		void add_site(int x, int neib[][4]); // adds a site to an existing bag
		void rem_site(int x, int neib[][4]); // removes a site from a bag
		void print_members ();
		double calc_det(int neib[][4]); // compute the bag determinant
		double get_trace_inv(int neib[][4]); // compute the trace  of the inverse bag Dirac operator
		bool in_bag(int x); // checks if a site belongs to a bag
		bool check_sites(int neib [][4]); //check if all sites in the bag are connected with at least 1 link
		bool check_overlap(bag b); //check if bag overlaps with bag b
		bag operator+(bag b); //add bag b to a bag
//		bag add_bags(bag a, bag b, int neib[][4]); //add bag a and bag b, the det gets calculated automatically
		void operator=(bag b); //assign bag b to a bag

};

//Default Constructor
bag::bag ()
{
	//bag_coords.push_back(0.0);
	det_D = 0.0;
}

//Default Destructor
bag::~bag ()
{

}

void bag::set_bag_coords(vector<int> vec)
{
	bag_coords = vec;
}

vector<int> bag::get_bag_coords()
{
	return bag_coords;
}

void bag::set_det_D(double det_M)
{
 	det_D = det_M;
}

double bag::get_det_D()
{
 return det_D;
}

// Add site x to the bag_coords
void bag::add_site(int x, int neib[][4])
{
//	bag_coords.push_back(x);
//	sort(bag_coords.begin(), bag_coords.end())

	if(bag_coords.empty())
	{
		bag_coords.push_back(x);
	}
	else
	{
		if(x > bag_coords.back())
		{
			bag_coords.push_back(x);
		}
		else
		{
			for(size_t i = 0; i < bag_coords.size(); i++)
			{
				if(x < bag_coords.at(i))
				{
					bag_coords.insert(bag_coords.begin()+i, x);
					break;
				}
			}
		}
	}
	set_det_D(calc_det(neib));
}

// Add site x to the bag_coords
void bag::rem_site(int x, int neib[][4])
{
	if(bag_coords.empty())
	{
		cout << "This bag is empty" << endl;
	}
	else{
	
		for(size_t i = 0; i < bag_coords.size(); i++)
		{
			if(bag_coords.at(i) == x)
			{
				bag_coords.erase(bag_coords.begin()+i);
				break;
			}
		}
	}
	set_det_D(calc_det(neib));
}

void bag::print_members()
{
  cout << endl << "This is a baryon bag!" << endl << endl;
  cout << "The coordinates of this bag are: ";
  for(auto it: bag_coords) cout << it << " ";
  cout << endl << endl;
  cout << "The fermion determinant is: " << det_D << endl;
 
}

double bag::calc_det(int neib[][4])
{
	vector<double> D(bag_coords.size()*bag_coords.size(), 0.0);
	
	for(size_t iy = 0; iy < bag_coords.size(); iy++)
	{	
		for(size_t ix = 0; ix < bag_coords.size(); ix++)
		{
			int k = ix + iy*bag_coords.size();
			int gamma = -1;
			if((bag_coords.at(ix) % Ns) % 2 == 0)  gamma = 1; // check if x1 of site is even or odd
		
//			cout << "--> " << k << endl;

			if(ix == iy) D.at(k) = 8.0*m*m*m;
			if(bag_coords.at(iy) == neib[bag_coords.at(ix)][0]) D.at(k) = gamma*1.0; //+ from kernel
			if(bag_coords.at(iy) == neib[bag_coords.at(ix)][1]) D.at(k) =  1.0; //+ from kernel, + from gamma
			if(bag_coords.at(ix) == neib[bag_coords.at(iy)][0]) D.at(k) = gamma*(-1.0); //- from kernel
			if(bag_coords.at(ix) == neib[bag_coords.at(iy)][1]) D.at(k) = -1.0; //- from kernel, + from gamma
			
//			cout << D.at(k) << endl;
		}
	}
	
	int s;
	gsl_matrix_view m = gsl_matrix_view_array (D.data(), bag_coords.size(), bag_coords.size());
	gsl_permutation * p = gsl_permutation_alloc (bag_coords.size());
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  
  cout << "det D = " << gsl_linalg_LU_det(&m.matrix, s) <<  endl;
  
  gsl_permutation_free(p);

// 	return gsl_linalg_LU_det(&m.matrix, s);
return 0;
}

double bag::get_trace_inv(int neib[][4])
{
	vector<double> D(bag_coords.size()*bag_coords.size(), 0.0);
	
	for(size_t iy = 0; iy < bag_coords.size(); iy++)
	{	
		for(size_t ix = 0; ix < bag_coords.size(); ix++)
		{
			int k = ix + iy*bag_coords.size();
			int gamma = -1;
			if((bag_coords.at(ix) % Ns) % 2 == 0)  gamma = 1; // check if x1 of site is even or odd
		
//			cout << "--> " << k << endl;

			if(ix == iy) D.at(k) = 8.0*m*m*m;
			if(bag_coords.at(iy) == neib[bag_coords.at(ix)][0]) D.at(k) = gamma*1.0; //+ from kernel
			if(bag_coords.at(iy) == neib[bag_coords.at(ix)][1]) D.at(k) =  1.0; //+ from kernel, + from gamma
			if(bag_coords.at(ix) == neib[bag_coords.at(iy)][0]) D.at(k) = gamma*(-1.0); //- from kernel
			if(bag_coords.at(ix) == neib[bag_coords.at(iy)][1]) D.at(k) = -1.0; //- from kernel, + from gamma
			
//			cout << D.at(k) << endl;
		}
	}
	
	int s;
	gsl_matrix_view m = gsl_matrix_view_array (D.data(), bag_coords.size(), bag_coords.size());
	gsl_matrix * m_inv = gsl_matrix_alloc(bag_coords.size(), bag_coords.size());
	gsl_permutation * p = gsl_permutation_alloc (bag_coords.size());
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_LU_invert(&m.matrix, p, m_inv);
  
  double tmp = 0.0;
  
  for(size_t i = 0; i < bag_coords.size(); i++)
  {
  	tmp += gsl_matrix_get(m_inv,i,i);
  }
  
  gsl_permutation_free(p);
  gsl_matrix_free(m_inv);
  
  return tmp;
}

bool bag::in_bag(int x)
{
	bool tmp = false;
	for(auto it: bag_coords)
	{
		if(it == x)
		{
			tmp = true;
			break;
		}
	}
	
	return tmp;
}

bool bag::check_sites(int neib [][4])
{
	bool bool_tmp = false;
	
	if(bag_coords.size() == 2)
	{
		int tmp = abs(bag_coords.at(0) - bag_coords.at(1));
		
		if(tmp == 1 || tmp == 4) bool_tmp = true;
		
		if(!bool_tmp) cout << "Alert: disconnected bag!" << endl;
	}
	else if(bag_coords.size() > 1)
	{
		for(auto it1: bag_coords)
		{	
			int neib0 =  neib[it1][0], neib1 =  neib[it1][1], neib2 =  neib[it1][2], neib3 =  neib[it1][3];
			bool_tmp = false;
			
			for(auto it2: bag_coords)
			{
				if(it2 == neib0)
				{
					bool_tmp = true;
				}
				else if(it2 == neib1)
				{
					bool_tmp = true;
				}
				else if(it2 == neib2)
				{
					bool_tmp = true;
				}
				else if(it2 == neib3)
				{
					bool_tmp = true;
				}
				
			}
			
			if(!bool_tmp) cout << "Alert: disconnected bag! (" << it1 << ")" << endl;
		}
	}
	
	return bool_tmp;
}

bool bag::check_overlap(bag b)
{
	bool bool_tmp = false;
	
	for(auto it_a: get_bag_coords())
	{
		for(auto it_b: b.get_bag_coords())
		{
			if(it_a == it_b)
			{
				bool_tmp = true;
				cout << "Bags overlap @ ( " << it_a << " )" << endl;
			}
		}
	}
	
	return bool_tmp;
}

//This is just the addition of two bags, in the simulation there should be a function merge_bags(a,b) which also deletes one old bag
bag bag::operator+(bag b)
{
	bag tmp_bag;
	vector<int> tmp_coords = get_bag_coords();
	
	for(auto it: b.get_bag_coords()) tmp_coords.push_back(it);
	sort(tmp_coords.begin(), tmp_coords.end());
	
	tmp_bag.set_bag_coords(tmp_coords);
	
	//calc det_D
	tmp_bag.set_det_D(0.0);
	
	return tmp_bag;
}

bag add_bags(bag a, bag b, int neib[][4])
{
	bag tmp_bag;
	vector<int> tmp_coords = a.get_bag_coords();
	
	for(auto it: b.get_bag_coords()) tmp_coords.push_back(it);
	sort(tmp_coords.begin(), tmp_coords.end());
	
	tmp_bag.set_bag_coords(tmp_coords);
	
	//calc det_D
	tmp_bag.set_det_D(tmp_bag.calc_det(neib));
	
	return tmp_bag;
}

void bag::operator=(bag b)
{
	set_bag_coords(b.get_bag_coords());
	set_det_D(b.get_det_D());
}

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


int main(){

/*
Goal: We want to have

vector<bag> bag_vec;

i = 0, ..., Nbag - 1, Nbag = # bags should be dynamical

*/

	bag a, b, d;
	int neib [V][4] = {0};
	
	neib_init(Nt, Ns, neib);
	
	//a.print_members();
	
	a.add_site(1, neib);
	b.add_site(2, neib);
	b.add_site(5, neib);
	
	a.check_overlap(b);
	
	bag c = a + b;
	
	//for(auto it: c.get_bag_coords()) cout << it << " ";	
		
	//(a+b).print_members();
	
	//cout << endl << endl;
	
//	a.set_bag_coords(vector<int>{1,5});
	
	c.calc_det(neib);
	
	cout << c.get_trace_inv(neib);
	
	bag e = add_bags(a, b, neib);
	
	e.print_members();
	
	e.calc_det(neib);
	
	cout << "+++++>" << e.get_det_D() << endl;
	cout << "+++++>" << e.get_trace_inv(neib) << endl;
	
	a.check_sites(neib);
	
//	a.rem_site(1);
//	
//	cout << "-> " << a.get_det_D() << endl;
//	
	

	return 0;
}



