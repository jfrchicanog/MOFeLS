//Algorithm which calculates aproximated well-spread Pareto fronts in MOCO using local search
//Version 1. Input maxt and TILIM
//Author: 
//Date of last modification : 25/10/2020



#include <list>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h> 
#include <algorithm>
#include <set>
#include <ctime>

//using namespace std;


#include "global_and_macros.h"
#include "moilp.h"
#include "heuristic2.h"



void Find_upper_bounds_that_dominate_z_and_create_Bj(MOILP &P, point &z, std::vector<std::multiset<Box *, Compare_value>> &XSet, std::list<Box *> &R, std::vector<std::list<Box *>> &Bj) {
	//Input : vector<multiset<Box *, Compare_value> XSet
	//Input : point z
	//Output : list<Box *> R  (boxes with upper bound dominating z. We use container to save the position (if neccesary))
	//Output : vector<list<Box *>> Bj	(Set Bj for all j=1,..,p. Used in RE algorithm)
	//Output : vector<multiset<Box *, Compare_value>> XSet (modified)   ( XSet = XSet - ( R U {Bj} ) )

	//Bj , j=1,..,p  is a list of boxes where the upper bound of its corresponding box have component j equals z.
	//The boxes from Bj are restored to XSet later in the algorithm
	
	int i, j;
	std::multiset<Box *>::iterator k;
	int coincidence_axis;
	int number_of_coincidences;
	bool dominates;
	int p = P.dimension;

	for (i = 0; i < XSet.size(); i++) { //For every tree in vector L
		k = XSet[i].begin();
		while (k != XSet[i].end()) {
			coincidence_axis = -1;
			number_of_coincidences = 0;
			dominates = true;
			for (j = 0; j < p; j++) {
				if (z[j] >(*k)->ub.at(j)) {
					dominates = false;
					j = p;
				}
				else if (z[j] == (*k)->ub.at(j)) {
					coincidence_axis = j;
					number_of_coincidences++;
				}
			}
			if (dominates) {
				if (coincidence_axis == -1) { //Strictly dominates
					R.push_back(*k);			//New box to R
					k = XSet[i].erase(k);
				}
				else if (number_of_coincidences == 1) {
					Bj[coincidence_axis].push_back(*k);	//New box to Bj
					k = XSet[i].erase(k);
				}
				else {
					++k;
				}
			}
			else {
				++k;
			}
		}
	}
#ifdef ON_SCREEN
	printf("|R|= %d", int(R.size()));
#endif
}

Box *PARTITION_BOX_full(MOILP &P, Box *B, int &coordinate, point &z, double *scaling) {
	int j;

	if (z[coordinate] == B->lb[coordinate])
		return NULL;

	int p = P.dimension;
	Box *newbox = new(Box);
	point lb(p), ub(p);

	for (j = 0; j < p; j++) {
		lb[j] = B->lb[j];
		ub[j] = B->ub[j];
	}
	ub[coordinate] = z[coordinate];
	newbox->set(lb, ub, coordinate, scaling, B);
		
	return (newbox);
}

Box *PARTITION_BOX_partition(MOILP &P, Box *B, int &coordinate, point &z, double *scaling) {
	
	point zz(P.dimension);
	int j;
	int p = P.dimension;

	Calculate_projection(z, zz, B);

	if (zz[coordinate] == B->lb[coordinate]) //Box is in a facet, then B_i(zz) is empty
		return NULL;

	Box *newbox = new(Box);
	point lb(p), ub(p);

	for (j = 0; j < p; j++) {
		if (j < coordinate) {
			lb[j] = B->lb[j];		ub[j] = B->ub[j];
		}
		else if (j == coordinate) {
			lb[j] = B->lb[j];		ub[j] = zz[j];
		}
		else {
			lb[j] = zz[j];		ub[j] = B->ub[j];
		}
	}
	newbox->set(lb, ub, coordinate, scaling, B);
	return (newbox);
}

void Split_boxes_and_create_Pj(MOILP &P, Input &input, point &z, std::list<Box *> &R, std::vector<std::list<Box *>> &Pj, double *scaling) {

	std::list<Box *>::iterator it;
	int j;
	int p = P.dimension;

	if (input.partition == "full") {
		for (it = R.begin(); it != R.end(); ++it) {
			for (j = 0; j < p; j++) {
				Box *u_j = NULL;
				u_j = PARTITION_BOX_full(P, *it, j, z, scaling);

				if (u_j != NULL) {
					Pj[j].push_back(u_j);
#ifdef ON_SCREEN
					show(*u_j, "new ");
#endif	
				}
			}
		}
	}
	else if (input.partition == "p-partition") {
		for (it = R.begin(); it != R.end(); ++it) {
			for (j = 0; j < p; j++) {
				Box *u_j = NULL;
				u_j = PARTITION_BOX_partition(P, *it, j, z, scaling);

				if (u_j != NULL) {
					Pj[j].push_back(u_j);
#ifdef ON_SCREEN
					show(*u_j, "new ");
#endif
				}
			}
		}
	}
}

void Filter_dominated_boxes_one_set(MOILP &P, Input &input, std::list<Box *> &Set, double *scaling) {
	std::list<Box *>::iterator it, it2;

	it = Set.begin();
	bool flag = false;
	int p = P.dimension;
		
	if (input.partition == "full") {
		while (it != Set.end()) {
			it2 = it;
			++it2;
			if (it2 == Set.end()) ++it;
			while (it2 != Set.end()) {
				flag = false;
				if (utils::vector_v1_lessorequal_v2((*it)->ub, (*it2)->ub, p)) { //if  (it <= it2)
#ifdef ON_SCREEN
					show(**it2, "modified");
#endif
					delete(*it);		P.stat.dominated_boxes++;
					it = Set.erase(it);
					it2 = Set.end();
					flag = true;
				}
				else if (utils::vector_v1_lessorequal_v2((*it2)->ub, (*it)->ub, p)) { //it2 <= it
#ifdef ON_SCREEN
					show(**it, "modified");
#endif
					delete(*it2);		P.stat.dominated_boxes++;
					it2 = Set.erase(it2);
				}
				else {
					++it2;
				}
			}
			if ((!flag) && (it != Set.end()))
				++it;
		}
	}
	else if (input.partition == "p-partition") {
		while (it != Set.end()) {
			it2 = it;
			++it2;
			if (it2 == Set.end()) ++it;
			while (it2 != Set.end()) {
				flag = false;
				if (utils::vector_v1_lessorequal_v2((*it)->ub, (*it2)->ub, p)) { //if  (it <= it2)
					Box::update_to_min_lb((*it2)->lb, (*it)->lb); //Updating it2->lb
					(*it2)->value = Box::get_value((*it2)->lb, (*it2)->ub, scaling);		//Recalculate the value when joining the two boxes
#ifdef ON_SCREEN
					show(**it2, "modified");
#endif
					delete(*it);		P.stat.dominated_boxes++;
					it = Set.erase(it);
					it2 = Set.end();
					flag = true;
				}
				else if (utils::vector_v1_lessorequal_v2((*it2)->ub, (*it)->ub, p)) { //it2 <= it
					Box::update_to_min_lb((*it)->lb, (*it2)->lb); //Updating it->lb
					(*it)->value = Box::get_value((*it)->lb, (*it)->ub, scaling);		//Recalculate the value when joining the two boxes
#ifdef ON_SCREEN
					show(**it, "modified");
#endif
					delete(*it2);		P.stat.dominated_boxes++;
					it2 = Set.erase(it2);
				}
				else {
					++it2;
				}
			}
			if ((!flag) && (it != Set.end()))
				++it;
		}
	}
}

void Filter_dominated_boxes_two_sets(MOILP &P, Input &input, std::list<Box *> &Set1, std::list<Box *> &Set2, double *scaling) {
	
	std::list<Box *>::iterator it, it2;

	it = Set1.begin();
	bool flag = false;
	int p = P.dimension;

	if (input.partition == "full") {
		while (it != Set1.end()) {
			it2 = Set2.begin();
			if (it2 == Set2.end()) ++it;
			while (it2 != Set2.end()) {
				flag = false;
				if (utils::vector_v1_lessorequal_v2((*it)->ub, (*it2)->ub, p)) { //it <= it2
#ifdef ON_SCREEN
					show(**it2, "modified");
#endif
					delete(*it);		P.stat.dominated_boxes++;
					it = Set1.erase(it);
					it2 = Set2.end();
					flag = true;
				}
				else if (utils::vector_v1_lessorequal_v2((*it2)->ub, (*it)->ub, p)) { //it2 <= it
#ifdef ON_SCREEN
					show(**it, "modified");
#endif
					delete(*it2);		P.stat.dominated_boxes++;
					it2 = Set2.erase(it2);
				}
				else {
					++it2;
				}
			}
			if ((!flag) && (it != Set1.end()))
				++it;
		}
	}
	else if (input.partition == "p-partition") {
		while (it != Set1.end()) {
			it2 = Set2.begin();
			if (it2 == Set2.end()) ++it;
			while (it2 != Set2.end()) {
				flag = false;
				if (utils::vector_v1_lessorequal_v2((*it)->ub, (*it2)->ub, p)) { //it <= it2
					Box::update_to_min_lb((*it2)->lb, (*it)->lb); //Updating it2->lb
					(*it2)->value = Box::get_value((*it2)->lb, (*it2)->ub, scaling);		//Recalculate the value when joining the two boxes
#ifdef ON_SCREEN
					show(**it2, "modified");
#endif
					delete(*it);		P.stat.dominated_boxes++;
					it = Set1.erase(it);
					it2 = Set2.end();
					flag = true;
				}
				else if (utils::vector_v1_lessorequal_v2((*it2)->ub, (*it)->ub, p)) { //it2 <= it
					Box::update_to_min_lb((*it)->lb, (*it2)->lb); //Updating it->lb
					(*it)->value = Box::get_value((*it)->lb, (*it)->ub, scaling);		//Recalculate the value when joining the two boxes
#ifdef ON_SCREEN
					show(**it, "modified");
#endif
					delete(*it2);		P.stat.dominated_boxes++;
					it2 = Set2.erase(it2);
				}
				else {
					++it2;
				}
			}
			if ((!flag) && (it != Set1.end()))
				++it;
		}
	}
}

void Filter_points_with_RE(MOILP &P, Input &input, std::vector<std::list<Box *>> &Pj, std::vector<std::list<Box *>> &Bj, double *scaling) {
	int j;

	for (j = 0; j < P.dimension; j++) { //Comparing Pj with elements of Pj 
		Filter_dominated_boxes_one_set(P, input, Pj[j], scaling);
	}
	for (j = 0; j < P.dimension; j++) { //Comparing Pj with elements of Bj 
		Filter_dominated_boxes_two_sets(P, input, Pj[j], Bj[j], scaling);
	}
}

void Insert_set_to_SET(MOILP &P, Input &input, std::vector<std::list<Box *>> &set, std::vector<std::multiset<Box *, Compare_value>> &XSET) {
	int j;
	std::list<Box *>::iterator it;
	
	if (input.set_of_boxes == "1") {
		for (j = 0; j < P.dimension; j++) {
			it = set[j].begin();
			while (it != set[j].end()) {
				XSET[0].insert(*it);
				it = set[j].erase(it);
			}
		}
	}
	else if (input.set_of_boxes == "alternate") {
		for (j = 0; j < P.dimension; j++) {
			it = set[j].begin();
			while (it != set[j].end()) {
				XSET[(*it)->direction].insert(*it);
				it = set[j].erase(it);
			}
		}
	}
}

void UPDATE_BOXES(MOILP &P, Input &input, point &z, std::vector<std::multiset<Box *, Compare_value>> &Set, std::list<Box *> &R, int &sizeL) {

	std::vector<std::list<Box *>> Bj(P.dimension);
	std::vector<std::list<Box *>> Pj(P.dimension);

	//First part. Calculate upper bounds that dominate z and Bj 
	Find_upper_bounds_that_dominate_z_and_create_Bj(P,z, Set, R, Bj);

	//Second part. Partition boxes and save them in Pj
	Split_boxes_and_create_Pj(P, input, z, R, Pj, P.APF.scaling);
	
	//Third part. Filtering solutions
	Filter_points_with_RE(P, input, Pj, Bj, P.APF.scaling);
	
	for (int i = 0; i < P.dimension; i++)		sizeL += int(Pj[i].size());
	sizeL -= int(R.size());

	//Fourth part. Update list
	Insert_set_to_SET(P,input, Pj, Set);
	Insert_set_to_SET(P,input, Bj, Set);
}

Box *Select_next_box(MOILP &P, Input &input, std::vector<std::multiset<Box *, Compare_value>> &L, int &counter) {
	
	if (input.set_of_boxes == "1") {
		return (*L[0].begin());
	}
	else if (input.set_of_boxes == "alternate") {
		//The next box to select varies in every component of L. We use a counter
		int index = counter % P.dimension;
		int maxcount = 0;
		while ( (L[index].size() == 0)  && (maxcount < P.dimension) ){	//In case of L(index) be empty
			counter++;
			index = counter % P.dimension;
			maxcount++;
		}
		if (maxcount == P.dimension) return NULL;
		counter++;
		return (*L[index].begin());
	}
	else {
		return NULL;
	}
}

void Add_new_APF_point(MOILP &P, double *x, point &z, double time) {

	//Insert a new solution, its image (x,z), and the time when the solution was obtained

	point xx(P.n_var), zz(P.dimension);
	solution A;

	for (int i = 0; i < P.n_var; i++)		xx.at(i) = x[i];
	for (int i = 0; i < z.size(); i++) zz.at(i) = z.at(i);
	A.x = xx;	A.z = zz;	A.time_point = time;
	
	P.APF.insert_without_domination(CPX_MIN, &A, P.dimension);

#ifdef ON_SCREEN
	show(&z, "z");		printf(" (total %d) ", int(P.APF.Set.size()));
#endif
}

void Delete_iterator_from_list(Input &input, Box *B, std::vector<std::multiset<Box *, Compare_value>> &L, int &sizeL) {
	if (input.set_of_boxes == "1") {
		L[0].erase(L[0].begin());
	}
	else if (input.set_of_boxes == "alternate") {
		L[B->direction].erase(L[B->direction].begin());
	}
	sizeL--;
}


bool improve(MOILP &P, int col_j, double x_j, char *sense, double &inc_c, double *slack , std::vector<double> &inc_slack ){

	double a_ij;	
	
	
	CPXgetobj(*P.cplex.env, *P.cplex.lp, &inc_c, col_j, col_j);
	if (inc_c == 0) 
		return false;		//No importa si la variable muta o no. No afecta en el coste
	else if (inc_c > 0) {
		if (round(x_j) == 0) 	
			return false;		//Introducir la variable en la solucion aumentaria el coste (estamos minimizando)
		else 					
			inc_c = -inc_c;		//En caso de retornar true, devolvemos el incremento en negativo, que sera la mejora en nuestra nueva solucion
	}
	else if ( (inc_c < 0) && (round(x_j) == 1) ) 
		return false;			//Quitar esa variable de la solucion empeora el objetivo, ya que aumentaria

	for(int i = 0 ; i < inc_slack.size() ; i++){		//Analizamos cada una de las restricciones
		inc_slack[i] = 0;
		CPXgetcoef(*P.cplex.env, *P.cplex.lp, i, col_j, &a_ij);

		switch (sense[i]){
			case 'L':
				if (round(x_j) == 0) {
					if ((a_ij >= 0) && (slack[i] <= a_ij))	
						return false;		//Meter la variable x_j hace que la restriccion no se cumpla	
					else
						inc_slack[i] = -a_ij;		//Si a_ij < 0 siempre puedo meter esa variable y actualizo slack. Si a_ij >= 0 pero slack[i] >= a_ij, también puedo meterla y el slack disminuye
				}
				else {
					if ( (a_ij < 0) && (slack[i] < abs(a_ij)) )
						return false;	//Quitar la variable aumenta el coste y sobrepaso el slack
					else
						inc_slack[i] = a_ij;	//Si a_ij >= 0 el slack aumenta, y si a_ij < 0 y slack_i >= abs(a_ij) puedo quitar esa variable de la solucion y el slack se ajusta
				}
				break;
			
			case 'G':
				if (round(x_j) == 0) {
					if ( (a_ij < 0) && (slack[i] > a_ij) )
						return false;		//Meter esa variable hace que el slack disminuya en mas de lo permitido. Restriccion infactible. Notese que slack <= 0 para restricciones de tipo >=
					else
 						inc_slack[i] = -a_ij; //Si a_ij >= 0, el slack mejora. Si a_ij < 0 y slack_i <= a_ij, puedo introducir la solucion y el slack se ajusta
				}
				else {
					if ( (a_ij >= 0) && (abs(slack[i])< a_ij) )
							return false;	
					else
							inc_slack[i] = a_ij;
				}
				break;
			
			case 'E':
				if (a_ij != 0) 
					return false;
			break;
		}
	}
	return true;
}




void print_vector(double *v, int n, std::string title){
	printf("\n%s " , title.c_str());
	for (int i = 0 ; i < n ; i++) printf("%.2f " , v[i]);
}


void Hill_climbing(MOILP &P, double *x, int m, int n, std::vector<int> &Index, double &obj, double *slack, char *sense){
	int i,j;
	double inc_c;
	std::vector<double> inc_slack(m);

	
	std::random_shuffle (Index.begin(), Index.end());

	for (j = 0 ; j < n ; j++){
		inc_c = 0;
		for (i = 0 ; i < m ; i++)	inc_slack[i] = 0;

		if (improve(P, Index[j], x[Index[j]], sense, inc_c, slack , inc_slack)){
				x[Index[j]] = 1 - x[Index[j]];		//Flip x_j
				obj += inc_c;						//Update objective
				for (i = 0 ; i < m ; i++){
					slack[i] += inc_slack[i];		//Update slack_j
				}
				j = -1 ;							//Restart
				//print_vector(x, P.n_var , "xx");
				//print_vector(slack, m, "slack");
				//printf(" ");
		}
		
	}
}

void Delete_boxes_in_R(std::list<Box *> &R){
	std::list<Box *>::iterator it = R.begin();
			while (it != R.end()) {
				delete(*it);
				it = R.erase(it);
			}
}

void RUN_ALGORITHM(MOILP &P, Input &input, std::vector<std::multiset<Box *, Compare_value>> &L, TIEMPO &t_ref) {
	
	std::srand ( unsigned ( std::time(0) ) );
	bool intime = true;						//Boolean to control time limit and max cplex iterations
		

	printf("\nExecuting problem %s with <%s %s %s %s> \n<TILIM %.1f>\n", input.name.c_str(), input.partition.c_str(), input.model.c_str(), 
		input.value.c_str(), input.set_of_boxes.c_str(), input.TILIM);
	printf("\nCalculating Non-dominated set...");

	Box *B0 = new (Box);
	B0->set(P.APF.Bound_I, P.APF.Bound_N, 0, P.APF.scaling, B0);		//Create initial box

	if (CPXsetintparam(*P.cplex.env, CPX_PARAM_MIPEMPHASIS, 1)) { warning = true;	printf("\nImpossible to set parameter MIPEMPHASIS to 1 in CPLEX");} //Fixed parameter	
	if (CPXsetintparam(*P.cplex.env, CPX_PARAM_INTSOLLIM, 1)) { warning = true;	printf("\nImpossible to set parameter INTSOLLIM to 1 in CPLEX");} //Fixed parameter
	if (CPXsetdblparam(*P.cplex.env, CPX_PARAM_TILIM, P.cplex.tilim)) { warning = true; 	printf("\nImpossible to set parameter TILIM in CPLEX"); }
	
	std::vector<int> Index(P.n_var);
	for (int i = 0 ; i < P.n_var ; i++) Index[i] = i;

	int m_ = P.n_const + P.dimension;
	double *x = (double *) malloc( P.n_var  * sizeof(double));
	double *slack = (double *) malloc( m_  * sizeof(double));
	char *sense = (char *) malloc ( m_ * sizeof(char));
	CPXgetsense(*P.cplex.env, *P.cplex.lp, sense, 0 , m_ -1);
	
	L[0].insert(B0);
	int sizeL = int(get_number_of_boxes_in_L(L));
	intime = Check_time_and_iterationslimit(t_ref, P.cplex);
	
	while ( (sizeL > 0) && (intime)) {
		int stat, count = 0;
		double obj;
		Box *B = NULL;
		point z(P.dimension);
		std::list<Box *> R;
	
		B = Select_next_box(P, input, L, count);

		#ifdef ON_SCREEN
			show(*B, "next box");
		#endif

		P.set_constraint_extreme_box(B->ub);
		P.cplex.set_time_for_solver(t_ref);
		
		
		obj = P.cplex.ILP_Solve_z(stat, &x, P.n_var);
		CPXgetslack(*P.cplex.env, *P.cplex.lp, slack, 0 , m_-1);
		//print_vector(slack, m_, "slack");
		//printf("  ");
		
		//if (stat == 104) || (stat == 107) || (stat == 101) || (stat == 102)
		if ( (stat == CPXMIP_SOL_LIM) || (stat == CPXMIP_TIME_LIM_FEAS) || (stat == CPXMIP_OPTIMAL) || (stat == CPXMIP_OPTIMAL_TOL) ){
			//print_vector(x,P.n_var, "x");
			Hill_climbing(P, x, m_, P.n_var, Index, obj, slack, sense);
			//print_vector(x,P.n_var, "x_");
			//print_vector(slack, m_, "slack");
			P.calculate_image_of_x(x, z);
			Add_new_APF_point(P, x, z, t_ref.value());
			UPDATE_BOXES(P, input, z, L, R, sizeL);
			Delete_boxes_in_R(R);
		}
		//if (stat == 103) || (stat == 119) || (stat == 108) 
		else if ( (stat  == CPXMIP_INFEASIBLE) || (stat  == CPXMIP_INForUNBD) || (stat == CPXMIP_TIME_LIM_INFEAS) ){
			Delete_iterator_from_list(input, B, L, sizeL);
			P.stat.empty_boxes++;
			delete(B);
		}
		else {
			printf("\nUnrecognizable status for the solver");
			exit(33);	
		}
		intime = Check_time_and_iterationslimit(t_ref, P.cplex);
	}

	free(x) ;
	free(slack);
	free(sense);
	
	//Finish execution. Free allocated memory from the remaining boxes
	P.stat.remaining_boxes = int(get_number_of_boxes_in_L(L));
	P.stat.dominated_points = int(P.APF.Discarded.size());
	destroy_list_boxes(L);
}



void Set_model(MOILP &P, Input &input, int argc, char **argv) {

	P.set_model(input.path, input.name, input.tmax, input.maxcplexiter, PARALLEL, THREADS, input.TILIM);

	

}

void Run_model(MOILP &P, Input &input) {
	TIEMPO t;
	std::vector<std::multiset<Box *, Compare_value>> L;				//In vector L we save all the boxes to analyze

	t.init();

	P.estimate_bounds_by_linear_relaxation(t, input.value);

	if (input.model == "chalmet") P.create_chalmet_model();
	else if (input.model == "tchebycheff") P.create_tchebycheff_model();

	if (input.set_of_boxes == "1") L.resize(1);
	else if (input.set_of_boxes == "alternate") L.resize(P.dimension);

	RUN_ALGORITHM(P, input, L, t);

	t.acum();
	P.stat.total_time = t.value();
	printf("DONE.\nEnd of execution after %lf seconds\n", P.stat.total_time);

	P.transform_pareto_front_if_neccesary();	//If original model was MAX, transform the solution
}

void Write_Results_file(MOILP &P, Input &input, std::string &date, std::string &hour, Metrics &M, bool flag) {
	FILE *fp;
	const char *name = "MOFeLS_results.txt";

	fp = fopen(name, "a+");
	char c;
	fread(&c, sizeof(c), 1, fp);	//Read the first character to check if the file is empty

	if (feof(fp)) { //If the file is empty, create the header
		fprintf(fp, "Date\t");
		fprintf(fp, "Hour\t");
		fprintf(fp, "Instance_name\t");
		fprintf(fp, "Maxt\t");
		fprintf(fp, "Algorithm\t");
		fprintf(fp, "|S|\t");
		fprintf(fp, "Iter\t"); 
		fprintf(fp, "Total_time\t");
		fprintf(fp, "TILIM\t"); 
		fprintf(fp, "Correct\t");
		fprintf(fp, "ONVGR\t");
		fprintf(fp, "GD\t");
		fprintf(fp, "IGD\t");
		fprintf(fp, "HV\t");
		fprintf(fp, "\n");
	}
	rewind(fp); //Rewind file

	//Write results in file
	fprintf(fp, "%s\t" ,date.c_str());
	fprintf(fp, "%s\t", hour.c_str());
	fprintf(fp, "%s\t", input.name.c_str());
	double tmax_ = input.tmax == MAX_DOUBLE ? 0 : input.tmax; fprintf(fp, "%.1f\t", tmax_);
	std::string algorithm = input.partition + "_" + input.model + "_" + input.value + "_" + input.set_of_boxes; fprintf(fp, "%s\t", algorithm.c_str());
	fprintf(fp, "%d\t", int(P.APF.Set.size()));
	fprintf(fp, "%d\t", int(P.cplex.n_iterations));
	fprintf(fp, "%.6f\t", P.stat.total_time);
	fprintf(fp, "%.3f\t", input.TILIM);
	warning == false ? fprintf(fp, "Yes\t") : fprintf(fp, "No\t");
	if (flag) {
		fprintf(fp, "%.3f\t", M.ONVGR);
		fprintf(fp, "%.10f\t", M.GD);
		fprintf(fp, "%.10f\t", M.IGD);
		fprintf(fp, "%.1f\t", M.HV);
	}
	fprintf(fp, "\n");
	fclose(fp);
}

bool Get_metrics(int argc, char **argv, MOILP &P, Input &input, Metrics &M) {

	std::string last_parameter = argv[argc-1];
	
	if ((last_parameter == "ONVG") || (last_parameter == "ONVGR") || (last_parameter == "GD") || (last_parameter == "IGD") || (last_parameter == "HV") ) {
		double target_irace;
		std::string pf_archive;
		bool wehavePF;
		pf_archive = input.path + ".pf";
		M.get_set(P.APF.Set);
		wehavePF = M.get_pf(pf_archive, P.dimension);
		
		if (wehavePF) {
			M.get_ONVG();
			M.get_ONVGR();
			M.get_GD();
			M.get_IGD();
#ifdef HV_FONSECA
			//M.get_HVf(DELTA_NADIR, P.original_sense, P.dimension);
#endif
#ifndef HV_FONSECA
			//M.get_HV(DELTA_NADIR, P.original_sense, P.dimension);		//M.get_HV(ref_point, type, p);		
#endif
			
			//Print selected metric on screen to use in IRACE
			if (last_parameter == "ONVG") target_irace = double(-M.ONVG);
			else if (last_parameter == "ONVGR") target_irace = -M.ONVGR;
			else if (last_parameter == "GD") target_irace = M.GD;
			else if (last_parameter == "IGD") target_irace = M.IGD;
			else if (last_parameter == "HV") target_irace = -M.HV;
			printf("\nOutput : %f", target_irace);
			return true;
		}
		else {
			return false;
		}	
	}
	return false;
}

void Print_output(MOILP &P, Input &input, int argc, char **argv, std::string folder) {
	
	std::string file_name, file1, file2, file3, date, hour;
	Metrics M;
	bool flag;

	utils::calculate_date_and_hour(date, hour);

		
#ifdef _WIN64
#include<direct.h>
	_mkdir(folder.c_str()); //Create a new folder. Return 0 if it does not exist, -1 otherwise
#elif __linux__
	mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

	file_name = folder;
	if (folder != "") file_name += "/";



	file_name = file_name + input.name + "_" + date + hour + std::to_string(P.stat.total_time);

	input.name += "_" + date + hour + std::to_string(P.stat.total_time);
	
	file1 = file_name;
	file2 = file_name + ".time";
	file3 = file_name + ".stat";
		
	P.APF.copy_set_to_file(file1, P.dimension);				//Print set of non-dominated points in a file
	P.APF.copy_set_and_time_to_file(file2, P.dimension);	//Print set of non-dominated with time in a file
	P.copy_stats_to_file(file3, input);						//Print statistics in a file

	flag = P.APF.check(P.dimension);
	if (!flag) {
		printf("\nTHE SET HAS DOMINATED POINTS");
		warning = true;
	}
	
	flag = Get_metrics(argc, argv, P, input, M);

	Write_Results_file(P, input, date, hour, M, flag);
}

int main(int argc, char **argv) {
//int main() {
	//main archive % tmax % maxcplexiter % partition % model % value % set_of_boxes % TILIM % [METRIC]
		
	Input input;								//Input parameters
	MOILP P;									//New Multiobjective object Model 

	Input_control(argc, argv, &input);			//Controlling that the input parameters are correct
		
	Set_model(P, input, argc,  argv);			//Set CPLEX parameters
	Run_model(P, input);						//Execute the algorithm
		
	Print_output(P, input, argc, argv, "PF");	//Print output into a folder
		
	return 0;
}
