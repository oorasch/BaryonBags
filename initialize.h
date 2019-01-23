#ifndef FILE_INITIALIZE
#define FILE_INITIALIZE
#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <vector>
#include "parameters.h"

//const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//std::ranlux48 generator(seed);

/*
 *
 * name: neib_init
 * @param N1, N2
 * @return filled neighbour field
 *
 */
void neib_init (const int &N1, const int &N2, int neib[][4]);

void init_md_config(const int neib[][4], std::vector<int>& site_type, std::vector<std::vector<int>>& dimer, std::vector<int>& bag_vec);

void config_init(std::vector<int>& site_type, std::vector<std::vector<int>>& dimer, std::vector<int>& bag_vec, std::vector<std::vector<bool>>& bag_site_vec);

/*void sig_init (std::vector<std::vector<int>>& sig_link);*/


#endif
