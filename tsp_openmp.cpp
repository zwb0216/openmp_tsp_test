#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <algorithm>
#include <climits>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <omp.h>



using namespace std;

int DEBUG = 0;
map<string, double> calculated_data;

	void print_map(map<string, vector<int> > mapa);
	void print_path(vector<string> path);
	void print_path_and_values(vector<string> path, map<string, vector<int> > mapa);
	void print_array_of_paths(vector<vector<string>> path_array);
	void print_array_of_paths_and_values(vector<vector<string>> path_array, map<string, vector<int> > mapa);

	double get_distance(string s1, string s2, map<string, vector<int> > mapa);
	double get_path_length(vector<string> nodes, map<string, vector<int> > mapa);
	double get_path_length(vector<string> nodes, map<string, vector<int> > mapa);
	double find_calculated_data(string id_1, string id_2);
	void add_to_calculated_data(string id_1, string id_2, double value);

	bool contains_cycles(vector<string> path);

	map<string, vector<int>> build_data_map();
	map<string, vector<int> > build_bigger_data_map();
	map<string, vector<int> > build_biggest_data_map();
	vector<string> get_nodes(map<string, vector<int> > mapa);
	vector<vector<string>> create_initial_population(map<string, vector<int>> mapa, unsigned int population_size);
	
	vector<string> shuffle(vector<string> nodes);
	vector<string> find_best_path(vector<vector<string>> population, map<string, vector<int> > mapa, unsigned int worst_possible_path_lenght);
	vector<vector<string>> pick_better_half(vector<vector<string>> population, map<string, vector<int> > mapa, unsigned int worst_possible_path_lenght);
	int get_random_num(int mod);
	vector<vector<string>> pick_two(vector<vector<string>> population);
	vector<vector<string>> crossover(vector<vector<string>> parents, map<string, vector<int> > mapa);
	vector<vector<string>> mutation(vector<vector<string>> children, map<string, vector<int> > mapa);
	void solve(vector<vector<string>> population_1, map<string, vector<int> > mapa, int iterations);

int main(){

	unsigned int population_size = 30;
	double mutation_threshold = 0.021;
	unsigned int iteration_threshold = 100;
	unsigned int worst_possible_path_lenght = UINT_MAX;
	srand(time(NULL));
	

	map<string, vector<int> > mapa = build_data_map();
	map<string, vector<int> > mapa2 = build_bigger_data_map();
	map<string, vector<int> > mapa3 = build_biggest_data_map();



	population_size = mapa2.size();

	vector<vector<string>> population_1 = create_initial_population(mapa3, population_size);

	clock_t begin_pt = clock();
	solve(population_1, mapa3, iteration_threshold);
	std::cout << "Time spent solving " << double(clock() - begin_pt) / CLOCKS_PER_SEC << endl;

	system("PAUSE");
	return 0;
	
}

void print_map(map<string, vector<int> > mapa){
	vector<string> v;
	for(map<string,vector<int>>::iterator it = mapa.begin(); it != mapa.end(); ++it) {
	  v.push_back(it->first);
	}
	for (unsigned int i = 0; i < v.size(); ++i){
		cout << v.at(i) << ": [ x: " << mapa[v.at(i)][0] << ", y: " << mapa[v.at(i)][1] << " ]"<< endl;
	}
}

void solve(vector<vector<string>> population_1, map<string, vector<int> > mapa, int iterations){

	double best_path_length = UINT_MAX;
	vector<string> best_path;
	unsigned int best_path_reached_in_iteration = 0;
	unsigned int iteration_counter = 0;
	unsigned int total_iterations = iterations;
	vector<vector<string>> local_population = population_1;


	while ( iterations > 0 ){
		vector<vector<string>> next_population;
		vector<vector<string>> better_half = pick_better_half(local_population, mapa, UINT_MAX);


		while (next_population.size() < population_1.size()){


			vector<vector<string>> parents = pick_two(better_half);


			vector<vector<string>> children = crossover(parents, mapa);

			children = mutation(children, mapa);

			bool cycles_in_daughter = contains_cycles(children.at(0));
			bool cycles_in_son = contains_cycles(children.at(1));


			if (children.at(0) != children.at(1)){
				if ( next_population.size() < 1 ){

					if (!cycles_in_daughter) next_population.push_back(children.at(0));
					if (!cycles_in_son) next_population.push_back(children.at(1));
				}
				else {

					int daughter_check = 0;
					int son_check = 0;
					for (unsigned int i = 0;  i < next_population.size(); ++i){
						if (children.at(0) == next_population.at(i)) ++daughter_check;
						if (children.at(1) == next_population.at(i)) ++son_check;
					}
					if (daughter_check < 1 && !cycles_in_daughter) next_population.push_back(children.at(0));
					if (son_check < 1 && !cycles_in_son) next_population.push_back(children.at(1));
				}
			}
			else {
				int daughter_check = 0;
				for (unsigned int i = 0; i < next_population.size(); ++i){
					if (children.at(0) == next_population.at(i)){
						++daughter_check;
					}
				}
				if (daughter_check < 1 && !cycles_in_daughter) next_population.push_back(children.at(0));
			}
		}

		vector<string> local_best_path = find_best_path(next_population, mapa, UINT_MAX);
		double local_best_path_length = get_path_length(local_best_path, mapa);

		cout << "Iteration £º" << iteration_counter+1 << ". Best path length: " << local_best_path_length << endl;

		if (local_best_path_length < best_path_length) {
			best_path_length = local_best_path_length;
			best_path = local_best_path;
			best_path_reached_in_iteration = total_iterations - iterations;
		}
		local_population = next_population;
		++iteration_counter;
		--iterations;
	}
	cout << "[FINAL] Best path length: " << best_path_length << " reached in " << best_path_reached_in_iteration << " iteration." << endl;

}

vector<vector<string>> mutation(vector<vector<string>> children, map<string, vector<int> > mapa){

	if (get_random_num(1000) > 21){
		return children;
	} 


	vector<vector<string>> mutated_children;
	vector<string> daughter = children.at(0);
	vector<string> son = children.at(1);

	int mutation_point_1 = get_random_num( (children.at(0)).size() );
	int mutation_point_2 = get_random_num( (children.at(0)).size() );
	while (mutation_point_1 == mutation_point_2){
		mutation_point_2 = get_random_num( (children.at(0)).size() );
	}

	unsigned int smaller = (mutation_point_1 < mutation_point_2)?(mutation_point_1):(mutation_point_2);
	unsigned int bigger = (mutation_point_1 < mutation_point_2)?(mutation_point_2):(mutation_point_1);

	vector<string> daughter_mutation;
	vector<string> son_mutation;
	for (unsigned int i = smaller; i < bigger; ++i){
		daughter_mutation.push_back(daughter.at(i));
		son_mutation.push_back(son.at(i));
	}

	reverse(daughter_mutation.begin(), daughter_mutation.end());
	reverse(son_mutation.begin(), son_mutation.end());

	unsigned int index = 0;
	for (unsigned int i = smaller; i < bigger; ++i){
		daughter[i] = daughter_mutation.at(index);
		son[i] = son_mutation.at(index);
		++index;
	}


	mutated_children.push_back(daughter);
	mutated_children.push_back(son);

	return mutated_children;
}

vector<vector<string>> crossover(vector<vector<string>> parents, map<string, vector<int> > mapa){

	vector<vector<string>> children;

	if ( get_random_num(100) > 70 ) return parents;

	int crossing_point_1 = get_random_num( (parents.at(0)).size() );
	int crossing_point_2 = get_random_num( (parents.at(0)).size() );
	while ( crossing_point_1 == crossing_point_2){
		crossing_point_2 = get_random_num( (parents.at(0)).size() );
	}

	unsigned int smaller = (crossing_point_1 < crossing_point_2)?(crossing_point_1):(crossing_point_2);
	unsigned int bigger = (crossing_point_1 < crossing_point_2)?(crossing_point_2):(crossing_point_1);

	vector<string> daughter;
	vector<string> son;


	for (unsigned int i = 0; i < (parents.at(0)).size(); ++i ){
		daughter.push_back(" ");
		son.push_back(" ");
	}

	vector<string> mother_data;
	vector<string> father_data;
	for (unsigned int i = smaller; i < bigger; ++i){
		daughter[i] = (parents.at(1)).at(i);
		mother_data.push_back((parents.at(1)).at(i));

		son[i] = (parents.at(0)).at(i);
		father_data.push_back((parents.at(0)).at(i));
	}


	vector<vector<string>> rules;
	for (unsigned int i = 0; i < mother_data.size(); ++i){
		vector<string> rule_1;
		vector<string> rule_2;
		rule_1.push_back(mother_data.at(i));
		rule_1.push_back(father_data.at(i));
		rule_2.push_back(father_data.at(i));
		rule_2.push_back(mother_data.at(i));
		rules.push_back(rule_1);
		rules.push_back(rule_2);
	}


	for (unsigned int i = 0; i < smaller; ++i){
		for (unsigned int j = 0; j < rules.size(); ++j){
			if ( (rules.at(j)).at(0) == (parents.at(0)).at(i) ){
				daughter[i] = (rules.at(j)).at(1);
			}
			if ( (rules.at(j)).at(0) == (parents.at(1)).at(i) ){
				son[i] = (rules.at(j)).at(1);
			}
		}
	}


	for (unsigned int i = bigger; i < daughter.size(); ++i){
		for (unsigned int j = 0; j < rules.size(); ++j){
			if ( (rules.at(j)).at(0) == (parents.at(0)).at(i) ){
				daughter[i] = (rules.at(j)).at(1);
			}
			if ( (rules.at(j)).at(0) == (parents.at(1)).at(i) ){
				son[i] = (rules.at(j)).at(1);
			}
		}
	}


	vector<string> unused_mommy_genes;
	vector<string> unused_daddy_genes;
	for (unsigned int i = 0; i < daughter.size(); ++i){
		if (daughter.at(i) == " ") unused_mommy_genes.push_back( (parents.at(0)).at(i) ); 
	}
	for (unsigned int i = 0; i < son.size(); ++i){
		if (son.at(i) == " ") unused_daddy_genes.push_back( (parents.at(1)).at(i) ); 
	}


	unsigned int index = 0;
	for (unsigned int i = 0; i < daughter.size(); ++i){
		if (daughter.at(i) == " "){
			daughter[i] = unused_daddy_genes.at(index);
			++index;
		}
	}
	index = 0;
	for (unsigned int i = 0; i < son.size(); ++i){
		if (son.at(i) == " "){
			son[i] = unused_mommy_genes.at(index);
			++index;
		}
	}

	children.push_back(daughter);
	children.push_back(son);
	return children;
}

vector<vector<string>> pick_two(vector<vector<string>> population)
{

	vector<vector<string>> pair;
	
	while (pair.size() != 2){


		unsigned int random = get_random_num( population.size() );
		vector<string> choosen = population.at(random);


		if (pair.size() == 0){
			pair.push_back(choosen);
		}

		else {
			if (pair.at(0) != choosen) pair.push_back(choosen);
		}
	}
	return pair;
}

int get_random_num(int mod)
{
	return rand() % mod;
}

vector<vector<string>> pick_better_half(vector<vector<string>> population, map<string, vector<int> > mapa, unsigned int worst_possible_path_lenght)
{

	vector<vector<string>> better_half;
	vector<vector<string>> local_population = population;

	for (unsigned int i = 0; i < population.size()/2; ++i){

		vector<string> best_one = find_best_path(local_population, mapa, worst_possible_path_lenght);
		better_half.push_back(best_one);

		for (unsigned int j = 0; j < local_population.size(); ++j){
			if (best_one == local_population.at(j)){

				local_population.erase(
					std::remove(local_population.begin(), local_population.end(), best_one), 
					local_population.end()
				);
			}
		}
	}
	return better_half;
}

vector<string> find_best_path(vector<vector<string>> population, map<string, vector<int> > mapa, unsigned int worst_possible_path_lenght)
{

	vector <string> best_path;
	double best_path_length = worst_possible_path_lenght;


	for (unsigned int i = 0; i < population.size(); ++i){
		double length = get_path_length(population.at(i), mapa);
		if (length < best_path_length){

			best_path_length = length;
			best_path = population.at(i);
		}
	}
	return best_path;
}

bool contains_cycles(vector<string> path)
{

	for (unsigned int i = 0; i < path.size(); ++i){
		int check = 0;
		for (unsigned int j = 0; j < path.size(); ++j){
			if (path.at(i) == path.at(j)){
	
				++check;
				if (check > 1) return true;
			}
		}
	}
	return false;
}

vector<string> get_nodes(map<string, vector<int> > mapa)
{

	vector<string> result;
	for (map<string, vector<int> >::iterator it = mapa.begin(); it != mapa.end(); ++it){
		result.push_back(it->first);
	}
	return result;
}

double get_distance(string s1, string s2, map<string, vector<int> > mapa)
{

	int x1 = mapa[s1].at(0);
	int y1 = mapa[s1].at(1);
	int x2 = mapa[s2].at(0);
	int y2 = mapa[s2].at(1);
	return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}

void print_path(vector<string> path){

	for (unsigned int i = 0; i < path.size(); ++i){
		cout << path.at(i) << " ";
	}
	cout << endl;
}

void print_path_and_values(vector<string> path, map<string, vector<int> > mapa){

	for (unsigned int i = 0; i < path.size(); ++i){
		cout << path.at(i) << " ";
	}
	cout << "\t" << get_path_length(path, mapa);
	cout << endl;
}

void print_array_of_paths(vector<vector<string>> path_array){

	for (unsigned int i = 0; i < path_array.size(); ++i){
		print_path(path_array.at(i));
	}
}

void print_array_of_paths_and_values(vector<vector<string>> path_array, map<string, vector<int> > mapa){

	for (unsigned int i = 0; i < path_array.size(); ++i){
		print_path_and_values(path_array.at(i), mapa);
	}
}

double old_get_path_length(vector<string> nodes, map<string, vector<int> > mapa){

	
	double sum = 0;
	for (unsigned int i = 0; i < nodes.size()-1; ++i){
		double partial_sum = find_calculated_data(nodes.at(i), nodes.at(i+1));
		if (partial_sum != -1){
			sum += partial_sum;
		}
		else {
			partial_sum = get_distance(nodes.at(i), nodes.at(i + 1), mapa);
			add_to_calculated_data(nodes.at(i), nodes.at(i+1), partial_sum);
			sum += partial_sum;
		}
	}
	return sum;
}

double get_path_length(vector<string> nodes, map<string, vector<int> > mapa){

	
	double sum = 0;
	int tid, nthreads, i;
	omp_set_num_threads(4);
	#pragma omp parallel shared(sum, calculated_data, mapa) private(tid, i)
	{
		tid = omp_get_thread_num();
	
		#pragma omp for schedule(dynamic)
		for (i = 0; i < nodes.size()-1; ++i){
			double partial_sum = find_calculated_data(nodes.at(i), nodes.at(i+1));
			if (partial_sum != -1){
				sum += partial_sum;
			}
			else {
				partial_sum = get_distance(nodes.at(i), nodes.at(i + 1), mapa);
				add_to_calculated_data(nodes.at(i), nodes.at(i+1), partial_sum);
				sum += partial_sum;
			}
		}
	}
	return sum;
}

vector<string> shuffle(vector<string> nodes){

	random_shuffle(nodes.begin(), nodes.end());
	return nodes;
}

vector< vector<string> > create_initial_population(map<string, vector<int> > mapa, unsigned int population_size){

	vector<vector<string>> population;
	vector<string> adam = get_nodes(mapa);


	while (population.size() < population_size){


		vector<string> being = shuffle(adam);
		if (population.size() < 1){
			population.push_back(being);
		}
		else {
			unsigned int check = 0;
			for (unsigned int i = 0; i < population.size(); ++i){
				if (being == population.at(i)){
					++check;
				}
			}
			if (check == 0){
				population.push_back(being);
			}
		}
	}

	for (unsigned int i = 0; i < population.size(); ++i){
		if ( contains_cycles( population.at(i))){
			cout << "[ERR] F2: Cycle detected at " << i << endl;
		}
	}
	return population;
}

void add_to_calculated_data(string id_1, string id_2, double value){

	string id = id_1 + "-" + id_2;
	string di = id_2 + "-" + id_1;
	calculated_data[id] = value;
	calculated_data[di] = value;
}

double find_calculated_data(string id_1, string id_2){

	string id = id_1 + "-" + id_2;
	double result = calculated_data[id];
	if (result != 0) {
		return result;
	}
	else {
		return -1;
	}
}

map<string, vector<int> > build_data_map(){

	vector<int> point;
	map<string, vector<int> > mapa;

	point.push_back(0);
	point.push_back(13);
	mapa["A"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(26);
	mapa["B"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(27);
	mapa["C"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(39);
	mapa["D"] = point;
	point.clear();

	point.push_back(2);
	point.push_back(0);
	mapa["E"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(13);
	mapa["F"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(19);
	mapa["G"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(25);
	mapa["H"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(31);
	mapa["I"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(37);
	mapa["J"] = point;
	point.clear();

	return mapa;
}

map<string, vector<int> > build_bigger_data_map(){

	vector<int> point;
	map<string, vector<int> > mapa;

	point.push_back(0);
	point.push_back(13);
	mapa["A"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(26);
	mapa["B"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(27);
	mapa["C"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(39);
	mapa["D"] = point;
	point.clear();

	point.push_back(2);
	point.push_back(0);
	mapa["E"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(13);
	mapa["F"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(19);
	mapa["G"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(25);
	mapa["H"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(31);
	mapa["I"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(37);
	mapa["J"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(43);
	mapa["K"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(8);
	mapa["L"] = point;
	point.clear();

	point.push_back(8);
	point.push_back(0);
	mapa["M"] = point;
	point.clear();

	point.push_back(9);
	point.push_back(10);
	mapa["N"] = point;
	point.clear();

	point.push_back(10);
	point.push_back(10);
	mapa["O"] = point;
	point.clear();

	point.push_back(11);
	point.push_back(10);
	mapa["P"] = point;
	point.clear();

	point.push_back(12);
	point.push_back(10);
	mapa["Q"] = point;
	point.clear();

	point.push_back(12);
	point.push_back(5);
	mapa["R"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(13);
	mapa["S"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(19);
	mapa["T"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(25);
	mapa["U"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(31);
	mapa["V"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(37);
	mapa["W"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(43);
	mapa["X"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(8);
	mapa["Y"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(11);
	mapa["Z"] = point;
	point.clear();
	return mapa;
}

map<string, vector<int> > build_biggest_data_map(){

	vector<int> point;
	map<string, vector<int> > mapa;



	point.push_back(0);
	point.push_back(13);
	mapa["001"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(26);
	mapa["002"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(27);
	mapa["003"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(39);
	mapa["004"] = point;
	point.clear();

	point.push_back(2);
	point.push_back(0);
	mapa["005"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(13);
	mapa["006"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(19);
	mapa["007"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(25);
	mapa["008"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(31);
	mapa["009"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(37);
	mapa["010"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(43);
	mapa["011"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(8);
	mapa["012"] = point;
	point.clear();

	point.push_back(8);
	point.push_back(0);
	mapa["013"] = point;
	point.clear();

	point.push_back(9);
	point.push_back(10);
	mapa["014"] = point;
	point.clear();

	point.push_back(10);
	point.push_back(10);
	mapa["015"] = point;
	point.clear();

	point.push_back(11);
	point.push_back(10);
	mapa["016"] = point;
	point.clear();

	point.push_back(12);
	point.push_back(10);
	mapa["017"] = point;
	point.clear();

	point.push_back(12);
	point.push_back(5);
	mapa["018"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(13);
	mapa["019"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(19);
	mapa["020"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(25);
	mapa["021"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(31);
	mapa["022"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(37);
	mapa["023"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(43);
	mapa["024"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(8);
	mapa["025"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(11);
	mapa["026"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(13);
	mapa["027"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(15);
	mapa["028"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(17);
	mapa["029"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(19);
	mapa["030"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(21);
	mapa["031"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(23);
	mapa["032"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(25);
	mapa["033"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(27);
	mapa["034"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(29);
	mapa["035"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(31);
	mapa["036"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(33);
	mapa["037"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(35);
	mapa["038"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(37);
	mapa["039"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(39);
	mapa["040"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(41);
	mapa["041"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(42);
	mapa["042"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(44);
	mapa["043"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(45);
	mapa["044"] = point;
	point.clear();

	point.push_back(25);
	point.push_back(11);
	mapa["045"] = point;
	point.clear();

	point.push_back(25);
	point.push_back(15);
	mapa["046"] = point;
	point.clear();

	point.push_back(25);
	point.push_back(22);
	mapa["047"] = point;
	point.clear();

	point.push_back(25);
	point.push_back(23);
	mapa["048"] = point;
	point.clear();

	point.push_back(25);
	point.push_back(24);
	mapa["049"] = point;
	point.clear();

	point.push_back(25);
	point.push_back(26);
	mapa["050"] = point;
	point.clear();

	point.push_back(25);
	point.push_back(28);
	mapa["051"] = point;
	point.clear();

	point.push_back(25);
	point.push_back(29);
	mapa["052"] = point;
	point.clear();

	point.push_back(25);
	point.push_back(9);
	mapa["053"] = point;
	point.clear();

	point.push_back(28);
	point.push_back(16);
	mapa["054"] = point;
	point.clear();

	point.push_back(28);
	point.push_back(20);
	mapa["055"] = point;
	point.clear();

	point.push_back(28);
	point.push_back(28);
	mapa["056"] = point;
	point.clear();

	point.push_back(28);
	point.push_back(30);
	mapa["057"] = point;
	point.clear();

	point.push_back(28);
	point.push_back(34);
	mapa["058"] = point;
	point.clear();

	point.push_back(28);
	point.push_back(40);
	mapa["059"] = point;
	point.clear();

	point.push_back(28);
	point.push_back(43);
	mapa["061"] = point;
	point.clear();

	point.push_back(28);
	point.push_back(47);
	mapa["062"] = point;
	point.clear();

	point.push_back(32);
	point.push_back(26);
	mapa["063"] = point;
	point.clear();

	point.push_back(32);
	point.push_back(31);
	mapa["064"] = point;
	point.clear();

	point.push_back(33);
	point.push_back(15);
	mapa["065"] = point;
	point.clear();

	point.push_back(33);
	point.push_back(26);
	mapa["066"] = point;
	point.clear();

	point.push_back(33);
	point.push_back(29);
	mapa["067"] = point;
	point.clear();

	point.push_back(33);
	point.push_back(31);
	mapa["068"] = point;
	point.clear();

	point.push_back(34);
	point.push_back(15);
	mapa["069"] = point;
	point.clear();

	point.push_back(34);
	point.push_back(26);
	mapa["070"] = point;
	point.clear();

	point.push_back(34);
	point.push_back(29);
	mapa["071"] = point;
	point.clear();

	point.push_back(34);
	point.push_back(31);
	mapa["072"] = point;
	point.clear();

	point.push_back(34);
	point.push_back(38);
	mapa["073"] = point;
	point.clear();

	point.push_back(34);
	point.push_back(41);
	mapa["074"] = point;
	point.clear();

	point.push_back(34);
	point.push_back(5);
	mapa["075"] = point;
	point.clear();

	point.push_back(35);
	point.push_back(17);
	mapa["076"] = point;
	point.clear();

	point.push_back(35);
	point.push_back(31);
	mapa["077"] = point;
	point.clear();

	point.push_back(38);
	point.push_back(16);
	mapa["078"] = point;
	point.clear();

	point.push_back(38);
	point.push_back(20);
	mapa["079"] = point;
	point.clear();

	point.push_back(38);
	point.push_back(30);
	mapa["080"] = point;
	point.clear();

	point.push_back(38);
	point.push_back(34);
	mapa["081"] = point;
	point.clear();

	point.push_back(40);
	point.push_back(22);
	mapa["082"] = point;
	point.clear();

	point.push_back(41);
	point.push_back(23);
	mapa["083"] = point;
	point.clear();

	point.push_back(41);
	point.push_back(32);
	mapa["084"] = point;
	point.clear();

	point.push_back(41);
	point.push_back(34);
	mapa["085"] = point;
	point.clear();

	point.push_back(41);
	point.push_back(35);
	mapa["086"] = point;
	point.clear();

	point.push_back(41);
	point.push_back(36);
	mapa["087"] = point;
	point.clear();

	point.push_back(48);
	point.push_back(22);
	mapa["088"] = point;
	point.clear();

	point.push_back(48);
	point.push_back(27);
	mapa["089"] = point;
	point.clear();

	point.push_back(48);
	point.push_back(6);
	mapa["090"] = point;
	point.clear();

	point.push_back(51);
	point.push_back(45);
	mapa["091"] = point;
	point.clear();

	point.push_back(51);
	point.push_back(47);
	mapa["092"] = point;
	point.clear();

	point.push_back(56);
	point.push_back(25);
	mapa["093"] = point;
	point.clear();

	point.push_back(57);
	point.push_back(12);
	mapa["094"] = point;
	point.clear();

	point.push_back(57);
	point.push_back(25);
	mapa["095"] = point;
	point.clear();

	point.push_back(57);
	point.push_back(44);
	mapa["096"] = point;
	point.clear();

	point.push_back(61);
	point.push_back(45);
	mapa["097"] = point;
	point.clear();

	point.push_back(61);
	point.push_back(47);
	mapa["098"] = point;
	point.clear();

	point.push_back(63);
	point.push_back(6);
	mapa["099"] = point;
	point.clear();

	point.push_back(64);
	point.push_back(22);
	mapa["100"] = point;
	point.clear();

	point.push_back(71);
	point.push_back(11);
	mapa["101"] = point;
	point.clear();

	point.push_back(71);
	point.push_back(13);
	mapa["102"] = point;
	point.clear();

	point.push_back(71);
	point.push_back(16);
	mapa["103"] = point;
	point.clear();

	point.push_back(71);
	point.push_back(45);
	mapa["104"] = point;
	point.clear();

	point.push_back(71);
	point.push_back(47);
	mapa["105"] = point;
	point.clear();

	point.push_back(74);
	point.push_back(12);
	mapa["106"] = point;
	point.clear();

	point.push_back(74);
	point.push_back(16);
	mapa["107"] = point;
	point.clear();

	point.push_back(74);
	point.push_back(20);
	mapa["108"] = point;
	point.clear();

	point.push_back(74);
	point.push_back(24);
	mapa["109"] = point;
	point.clear();

	point.push_back(74);
	point.push_back(29);
	mapa["110"] = point;
	point.clear();

	point.push_back(74);
	point.push_back(35);
	mapa["111"] = point;
	point.clear();

	point.push_back(74);
	point.push_back(39);
	mapa["112"] = point;
	point.clear();

	point.push_back(74);
	point.push_back(6);
	mapa["113"] = point;
	point.clear();

	point.push_back(77);
	point.push_back(21);
	mapa["114"] = point;
	point.clear();

	point.push_back(78);
	point.push_back(10);
	mapa["115"] = point;
	point.clear();

	point.push_back(78);
	point.push_back(32);
	mapa["116"] = point;
	point.clear();

	point.push_back(78);
	point.push_back(35);
	mapa["117"] = point;
	point.clear();

	point.push_back(78);
	point.push_back(39);
	mapa["118"] = point;
	point.clear();

	point.push_back(79);
	point.push_back(10);
	mapa["119"] = point;
	point.clear();

	point.push_back(79);
	point.push_back(33);
	mapa["120"] = point;
	point.clear();

	point.push_back(79);
	point.push_back(37);
	mapa["121"] = point;
	point.clear();

	point.push_back(80);
	point.push_back(10);
	mapa["122"] = point;
	point.clear();

	point.push_back(80);
	point.push_back(41);
	mapa["123"] = point;
	point.clear();

	point.push_back(80);
	point.push_back(5);
	mapa["124"] = point;
	point.clear();

	point.push_back(81);
	point.push_back(17);
	mapa["125"] = point;
	point.clear();

	point.push_back(84);
	point.push_back(20);
	mapa["126"] = point;
	point.clear();

	point.push_back(84);
	point.push_back(24);
	mapa["127"] = point;
	point.clear();

	point.push_back(84);
	point.push_back(29);
	mapa["128"] = point;
	point.clear();

	point.push_back(84);
	point.push_back(34);
	mapa["129"] = point;
	point.clear();

	point.push_back(84);
	point.push_back(38);
	mapa["130"] = point;
	point.clear();

	point.push_back(84);
	point.push_back(6);
	mapa["131"] = point;
	point.clear();

	point.push_back(107);
	point.push_back(27);
	mapa["132"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(13);
	mapa["133"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(26);
	mapa["134"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(27);
	mapa["135"] = point;
	point.clear();

	point.push_back(0);
	point.push_back(39);
	mapa["136"] = point;
	point.clear();

	point.push_back(2);
	point.push_back(0);
	mapa["137"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(13);
	mapa["138"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(19);
	mapa["139"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(25);
	mapa["140"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(31);
	mapa["141"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(37);
	mapa["142"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(43);
	mapa["143"] = point;
	point.clear();

	point.push_back(5);
	point.push_back(8);
	mapa["144"] = point;
	point.clear();

	point.push_back(8);
	point.push_back(0);
	mapa["145"] = point;
	point.clear();

	point.push_back(9);
	point.push_back(10);
	mapa["146"] = point;
	point.clear();

	point.push_back(10);
	point.push_back(10);
	mapa["147"] = point;
	point.clear();

	point.push_back(11);
	point.push_back(10);
	mapa["148"] = point;
	point.clear();

	point.push_back(12);
	point.push_back(10);
	mapa["149"] = point;
	point.clear();

	point.push_back(12);
	point.push_back(5);
	mapa["150"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(13);
	mapa["151"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(19);
	mapa["152"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(25);
	mapa["153"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(31);
	mapa["154"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(37);
	mapa["155"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(43);
	mapa["156"] = point;
	point.clear();

	point.push_back(15);
	point.push_back(8);
	mapa["157"] = point;
	point.clear();

	point.push_back(18);
	point.push_back(11);
	mapa["158"] = point;
	point.clear();

	return mapa;
}
