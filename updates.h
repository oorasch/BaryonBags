#ifndef FILE_UPDATES
#define FILE_UPDATES

void do_updates(const int& sweeps, const int neib[][4], std::vector<int>& site_type, std::vector<std::vector<int>>& dimer, std::vector<int>& bag_vec, std::vector<std::vector<bool>>& bag_sites, std::vector<double>& bag_det);

void B_bar_update(int& site_type_x, int& site_type_y, int& d);

void B_update(const int& x, const int& y, const int neib [][4], int& site_type_x, int& site_type_y, int& d, std::vector<int>& bag_vec, std::vector<std::vector<bool>>& bag_sites, std::vector<double>& bag_det);

void check_local_bag_occupation(const int& site, const int neib[][4], const std::vector<int>& bag_vec, std::vector<std::vector<int>>& local_bags);

#endif
