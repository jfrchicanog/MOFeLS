#pragma once
//Author: 
//Date of last modification : 25/10/2020


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <stdbool.h>	
#include <math.h>		
#include <vector>		
#include <list>




#ifdef ON_SCREEN
void show(point *A) {
	printf("(");
	for (int i = 0; i < A->size() - 1; i++) {
		printf("%2.0f,", A->at(i));
	}
	printf("%2.0f) ", A->at(A->size() - 1));
}

void show(point *A, std::string title) {
	printf("\n%s = (", title.c_str());
	for (int i = 0; i < A->size() - 1; i++) {
		printf("%2.0f,", A->at(i));
	}
	printf("%2.0f) ", A->at(A->size() - 1));
}

void show(int stat) {
	if (stat == 101) printf("\nEXACT");
	else if ((stat == 102) || (stat == 104) || (stat == 107)) printf("\nAPPROX");
	else if (stat == 108) printf("TILIM INFEASIBLE");
	else if ((stat != 103) && (stat != 119)) printf("UNRECOGNIZABLE");
}
void show(Box &B, std::string title) {
	int p = int(B.ub.size());
	printf("\n%s = %d {", title.c_str(), B.id);
	show(&B.lb);
	show(&B.ub);
	printf("} %0.3lf", B.value);
	printf(" p= ");
	//if (B.box_id != 0) printf("%d", B.parent->box_id);
	if ((B.parent->id < 0) || (B.parent->id > box_id)) printf(" ");
	else printf("%d ", B.parent->id);
	printf("p_sib= %d", B.direction);
	//printf(" c= ");
	//for (int i = 0; i < p; i++) {
	//if (B.children.at(i) == NULL)
	//printf("X ");
	//else
	//printf("%d ", B.children.at(i)->box_id);
	//}
}
#endif


double max(const double &a1, const double &a2) {
	return (a1 > a2 ? a1 : a2);
}

bool Check_time_and_iterationslimit(TIEMPO &t, CPLEX &cplex) {
	//Return false if time limit or the desired number of iterations are reached. Return true, otherwise
	t.acum();
	return ((t.value() < cplex.t_max) && (cplex.n_iterations < cplex.max_iter));
}

int get_number_of_boxes_in_L(std::vector<std::multiset<Box *, Compare_value>> &L) {
	int val = 0;
	for (int i = 0; i < L.size(); i++) {
		val += int(L[i].size());
	}
	return (val);
}


void error_input(void) {
	printf("\nTo execute MOFeLS write ./MOFeLS (FILE) (MAXT) (MAXIT) (arg4) (arg5) (arg6) (arg7) (TILIM) [METRIC]");
	printf("\n(FILE) is the name of the model. There must be two existing archives, (FILE).obj and (FILE).lp");
	printf("\n(MAXT) is the maximum total execution time in seconds (0 for unlimited time).");
	printf("\n(MAXIT) is the maximum number of CPLEX iterations (0 for unlimited size).");
	printf("\n(arg4) is the type of partition. Type 'full' or 'p-partition'.");
	printf("\n(arg5) is the parameterization model. Type 'chalmet' or 'tchebycheff'.");
	printf("\n(arg6) is the box-value. Type 'volume' or 'scaled'.");
	printf("\n(arg7) is the number of set of boxes used. Type '1' or 'alternate'.");
	printf("\n(TILIM) is the maximum time (in seconds) used by CPLEX to find a feasible solution.");
	printf("\n\n[METRIC] is an optional parameter used as target in IRACE. Type 'ONVG', 'ONVGR', 'GD', 'IGD' or 'HV'. In case of a selected metric, all of them will be calculated. A complete Pareto Set must be provided in the same folder, with name (FILE).pf.\n\n");
	exit(EXIT_FAILURE);
}


void Input_control(int argc, char **argv, Input *input) {
	//Controlling the input parameters. In case of error, exit.
	FILE *fp;


	if (argc < 9) {
		error_input(); 
	}

	//Check _obj file and .lp file
	std::string str = argv[1];
	str += ".obj";
	fp = fopen(str.c_str(), "r");
	if (fp == NULL)		
	utils::error_open_file(str.c_str());
	else 
	fclose(fp);

	str = argv[1];
	str += ".lp";
	fp = fopen(str.c_str(), "r");
	if (fp == NULL)		utils::error_open_file(str.c_str());
	else fclose(fp);


	//Maximum execution time
	double tmax = atof(argv[2]);
	if (tmax < 0) error_input();


	//Limit size in CPLEX iterations
	int limitcplexiter = atoi(argv[3]);
	if (limitcplexiter < 0)		error_input();


	//Partition
	std::string type = argv[4];
	if ((type != "full") && (type != "p-partition"))  error_input();

	//Problem model
	std::string model = argv[5];
	if ((model != "chalmet") && (model != "tchebycheff"))	error_input();

	//Value
	std::string value = argv[6];
	if ((value != "volume") && (value != "scaled"))	error_input();

	//Set of boxes
	std::string set_b = argv[7];
	if ((set_b != "alternate") && (set_b != "1")) error_input();

	//TILIM
	double tilim = atof(argv[8]);
	if (tilim <= 0)	error_input();

	
	input->path = argv[1];
	input->name = utils::get_name_string(input->path);
	input->tmax = atof(argv[2]);
	input->maxcplexiter = atoi(argv[3]);
	input->partition = argv[4];
	input->model = argv[5];
	input->value = argv[6];
	input->set_of_boxes = argv[7];
	input->TILIM = atof(argv[8]);
}

void Transfer_boxes_from_X_to_L(Input &input, std::list<Box *> &X, std::vector<std::multiset<Box *, Compare_value>> &L) {
	std::list<Box *>::iterator it = X.begin();

	if (input.set_of_boxes == "1") {
		while (it != X.end()) {
			L[0].insert(*it);
			it = X.erase(it);
		}
	}
	else if (input.set_of_boxes == "alternate") {
		while (it != X.end()) {
			L[(*it)->direction].insert(*it);
			it = X.erase(it);
		}
	}
}

void Transfer_boxes_from_L_to_X(std::vector<std::multiset<Box *, Compare_value>> &L, std::list<Box *> &X) {
	for (int j = 0; j < L.size(); j++) {
		for (std::multiset<Box *, Compare_value>::iterator it = L[j].begin(); it != L[j].end(); ++it) {
			X.push_back(*it);
		}
		L[j].clear();
	}


}

void Transfer_boxes_from_X_to_Y_eliminating_domination(MOILP &P, Input &input, std::list<Box *> &X, std::list<Box *> &Y) {

	std::list<Box *> Discarded;

	std::list<Box *>::iterator it = X.begin();
	if (input.partition == "full") {
		while (it != X.end()) {
			(*it)->insert_without_domination(&Y, &Discarded, P.dimension);
			it = X.erase(it);
		}
	}
	else if (input.partition == "p-partition") {
		while (it != X.end()) {
			(*it)->insert_without_domination_and_join(&Y, &Discarded, P.dimension);
			it = X.erase(it);
		}
	}

	P.stat.dominated_boxes += (int)Discarded.size();
	it = Discarded.begin();
	while (it != Discarded.end()) {
		delete(*it);
		it = Discarded.erase(it);
	}

}

void destroy_list_boxes(std::vector<std::multiset<Box *, Compare_value>> &L) {
	std::multiset<Box *>::iterator it;
	for (int i = 0; i < L.size() - 1; i++) {
		it = L[i].begin();
		while (it != L[i].end()) {
			delete(*it);
			it = L[i].erase(it);
		}
	}
}

void destroy_list_boxes(std::list<Box *> &Set) {
	std::list<Box *>::iterator it2;

	it2 = Set.begin();
	while (it2 != Set.end()) {
		delete(*it2);
		it2 = Set.erase(it2);
	}
}

void Calculate_projection(point &z, point &zz, Box *B) {
	int i;
	int p = int(z.size());
	for (i = 0; i < p; i++) {
		if (z[i] < B->lb[i]) zz[i] = B->lb[i];
		else zz[i] = z[i];
	}
}
