#ifndef FILE_DEBUG
#define FILE_DEBUG
#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <vector>
#include "parameters.h"

//const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//std::ranlux48 generator(seed);

void print_md_occ(const std::vector<int>& site_type, const std::vector<std::vector<int>>& dimer);

void print_lattice_occ(const std::vector<int>& site_type, const std::vector<std::vector<int>>& dimer);

void check_occupation (const int neib[][4], const std::vector<int>& site_type, const std::vector<std::vector<int>>& dimer);

void fill_bag_vec(const int neib[][4], const std::vector<int>& monomer, const std::vector<std::vector<int>>& dimer, std::vector<std::vector<int>>& bag_vec);

bool is_neib(const int neib[][4], const std::vector<int>& bag_coords, const int& x);


#endif
