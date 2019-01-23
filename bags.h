#ifndef BAGS_H
#define BAGS_H

void make_bag(const int& site, std::vector<int>& bag_vec, std::vector<std::vector<bool>>& bag_sites, std::vector<double>& bag_det);

void del_bag_1(const int& site, std::vector<int>& bag_vec, std::vector<std::vector<bool>>& bag_sites, std::vector<double>& bag_det);

void add_1_site(const int& site, const int& bag_label, std::vector<bool>& bag_sites_bag, std::vector<double>& bag_det);

void rem_1_site(const int& site, const int& bag_label, std::vector<bool>& bag_sites, std::vector<double>& bag_det);

double calc_det(const int neib [][4], const std::vector<bool>& bag);

#endif
