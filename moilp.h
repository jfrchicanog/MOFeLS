#pragma once
//Author: 
//Date of last modification : 25/10/2020

//*******************************
//LIBRARIES
//*******************************

#include <ilcplex/cplex.h>
#include <string>
#include <list>
#include <vector>
#include <time.h>
#include <set>
#include <math.h>

//*******************************
//DECLARATIONS
//*******************************
class TIEMPO;
class utils;
class Box;
class Metrics;
class Compare_value;
class CPLEX;
class Solution;
class MOILP;



//*******************************
//TYPEDEFS
//*******************************
typedef  std::vector<double>  point;

//using namespace std;


//*******************************
//STRUCTURES
//*******************************

struct solution { //Save the solutions (x,z) where z = f(x) and the time when obtaining the solution
	point x;
	point z;
	double time_point;
};

struct Input {
	std::string path;
	std::string name;
	double tmax;
	int maxcplexiter;
	std::string partition, model, value, set_of_boxes;
	double TILIM;
};

struct Stat {
public:
	int empty_boxes = 0;
	int remaining_boxes = 0;
	double total_time = 0.0;
	int dominated_boxes = 0;
	int dominated_points = 0;
};




//*******************************
//CLASSESS
//*******************************

#ifdef _WIN64
	#include <direct.h>

class TIEMPO {
public:
	clock_t  t, t2;
	int tipo;
public:
	TIEMPO() {	tipo = 0;	} //Constructor
	void TIEMPO::init() {
		t = clock();		t2 = 0;
	}
	void TIEMPO::acum() { //Acummulated time 
		clock_t t3 = clock();
		t2 += clock() - t;		t = t3;
	}
	double TIEMPO::value() { //Transform into seconds
		return ((double)(t2 / double(CLOCKS_PER_SEC)));
	}
	int gettype() {	return tipo;	}
	void chgtype(int newt) {	tipo = newt;	}

	~TIEMPO() { } //Destructor
};

#elif __linux__
#include <sys/stat.h> //Para mkdir en Linux
class TIEMPO {
	timespec t1, t2;
	int tipo;
	double time;
private:
	timespec diff(timespec start, timespec end) {
		timespec temp;
		if ((end.tv_nsec - start.tv_nsec)<0) {
			temp.tv_sec = end.tv_sec - start.tv_sec - 1;
			temp.tv_nsec = 1000000000.0 + end.tv_nsec - start.tv_nsec;
		}
		else {
			temp.tv_sec = end.tv_sec - start.tv_sec;
			temp.tv_nsec = end.tv_nsec - start.tv_nsec;
		}
		return temp;
	}
	timespec add(timespec accum, timespec now) {
		timespec temp;
		if (accum.tv_nsec + now.tv_nsec >= 1E9) {
			temp.tv_sec = accum.tv_sec + now.tv_sec + 1;
			temp.tv_nsec = accum.tv_nsec + now.tv_nsec - 1E9;
		}
		else {
			temp.tv_sec = accum.tv_sec + now.tv_sec;
			temp.tv_nsec = accum.tv_nsec + now.tv_nsec;
		}
		return temp;
	}
public:
	TIEMPO() {	//Create object TIEMPO and initializes to one of the four types
		tipo = CLOCK_REALTIME;
		//tipo = CLOCK_MONOTONIC;
		//tipo = CLOCK_PROCESS_CPUTIME_ID;
		//tipo = CLOCK_THREAD_CPUTIME_ID;
	} //Constructor

	void init() { //Restart
		time = 0.0;
		clock_gettime(tipo, &t1);
		if (tipo == CLOCK_THREAD_CPUTIME_ID) {
			t1.tv_sec = 0.0;			t1.tv_nsec = 0.0;
		}
	}
	void acum() { //Acummulate 
		timespec aux;

		clock_gettime(tipo, &t2);

		if (tipo == CLOCK_THREAD_CPUTIME_ID) {
			aux = diff(t1, t2);			t1 = add(t1, aux);
			time += t1.tv_sec + t1.tv_nsec / 1E9;
		}
		else if (tipo == CLOCK_REALTIME) {
			time += (double)(1.0*(1.0*t2.tv_nsec - t1.tv_nsec*1.0)*1e-9 + 1.0*t2.tv_sec - 1.0*t1.tv_sec);
		}
		t1 = t2;
	}
	double value() { //Transform into seconds
		return (time);
	}
	int gettype() {	return tipo;}

	void chgtype(int newt) {
		//0 = CLOCK_REALTIME
		//1 = CLOCK_MONOTONIC
		//2 = CLOCK_PROCESS_CPUTIME_ID
		//3 = CLOCK_THREAD_CPUTIME_ID
		tipo = newt;
	}
	~TIEMPO() { } //Destructor
};
#elif __APPLE__
//APPLE
#elif __unix__
//Unix
#elif defined(_POSIX_VERSION)
// POSIX
#else
#   error "Unknown compiler"
#endif


class utils {
public:
	utils() {}

public:
	static std::string get_name_string(const std::string &name_) {
		//Take the name of an archive, ommiting the previous path of the archive
		//Example: Input   ./KP/data/hello.lp		Output	hello.lp

		std::size_t pos;
		std::string name = name_;

		pos = name.find("/");
		while (pos != -1) {
			name = name.substr(pos + 1);
			pos = name.find("/");
		}
		return name;
	}
	static void error_CPLEX_status(const int &code) {
		printf("\nCPLEX error status = %d", code);
		exit(EXIT_FAILURE);
	}
	static int vector_v1_lessorequal_v2(const std::vector<double> &v1, const std::vector<double> &v2, const int &p) {
		//Compare two p-dimensional vectors v1 and v2. Return 1 if v1 <= v2 and return 0, otherwise
		for (int index = 0; index < p; index++) {
			if (v1[index] > v2[index]) return 0;
		}
		return 1;
	}
	static void error_open_file(const char *name) {
		printf("\nImpossible to open file %s\n", name);
		exit(EXIT_FAILURE);
	}

	static void calculate_date_and_hour(std::string &fecha, std::string &hora) {
		struct tm *newtime;
		time_t long_time;
		time(&long_time);
		newtime = localtime(&long_time);

		if (newtime->tm_mday < 10) 		fecha += ("0");
		fecha += std::to_string(newtime->tm_mday);
		if (newtime->tm_mon < 9)		fecha += ("0");
		fecha += std::to_string(newtime->tm_mon + 1);		//Month
		fecha += std::to_string(newtime->tm_year + 1900);	//Year

		if (newtime->tm_hour < 10) hora += ("0");
		hora += std::to_string(newtime->tm_hour);			//Hours
		if (newtime->tm_min < 10) hora += ("0");
		hora += std::to_string(newtime->tm_min);			//Minutes
		if (newtime->tm_sec < 10) hora += ("0");
		hora += std::to_string(newtime->tm_sec);			//Seconds
	}

	static int compare_lex(point *v1, point *v2, int index, int &p) {
	//Compare two points v1 and v2 in lexicographic order.
	//Return -2 if v1 <_lex v2 and not v1 < v2
	//Return -1 if v1 < v2 
	//Return 0 if v1 = v2
	//Return 1  if v1 > v2 
	//Return 2 if v1 >_lex v2 and not v1 > v2

		if (v1->at(index) < v2->at(index)) {
		for (int i = index + 1; i < p; i++) {
			if (v1->at(i) > v2->at(i))
				return -2;
		}
		return -1;
		}
		else if (v1->at(index) > v2->at(index)) {
		for (int i = index + 1; i < p; i++) {
			if (v1->at(i) < v2->at(i))
				return 2;
		}
		return 1;
		}
		else {
		if (index < p - 1)
			return (compare_lex(v1, v2, index + 1, p));
		else
			return 0; //Equal vectors
		}
	}


	static int insert_without_domination(int type, point *v, std::list<point> &Set, std::list<point> &Discarded, int &p) {
		//Insert a new nondominated point z in lexicographic order into the list L.
		//If type = 1, use <_lex		//If type = -1, use >_lex
		//Return 0 if success
		//Return -1 if the solution is dominated (thus, not included in the list)
		//Return a >= 1 , where a is the number of solutions which are deleted from the list, because they are dominated by the new solution

		//Outputs of compare_lex. //Compare two p-dimensional vectors v1 and v2 in lexicographic order.
		//Return -2 if v1 < v2; Return -1 if v1 <_lex v2 and not v1 < v2; Return 0 if v1 = v2; Return 1  if v1 > v2; Return 2 if v1 >_lex v2 and not v1 > v2

		std::list<point>::iterator it1, it2, aux;
		int n_eliminated = 0;
		bool finish = false;

		//Introduce the new point in the front of the list. Then compare it with every element on the list, one by one until its lexicographical position is found.
		Set.push_front(*v);
		it1 = it2 = Set.begin();

		if (Set.size() == 2) { //There is only one element to compare with
			std::advance(it2, 1);
			int aa = utils::compare_lex(&(*it1), &(*it2), 0, p);

			if (type == CPX_MIN) {
				if (aa == 1) { //it1 = v >  it2
					Discarded.push_back(Set.front());
					Set.pop_front();
					return -1;
				}
				else if (aa == 2) { //it1 = v >_lex it2
					Set.push_back(*it1);				Set.pop_front();
				}
				else if (aa == 0) { //The two vectors are equal   it1 = it2
					Discarded.push_back(Set.front());
					Set.pop_front();
					return -1;
				}
				else if (aa == -2) {//it1 = v <_lex it2. They are ordered

				}
				else if (aa == -1) { //it1 <  it2. We eliminate it2
					Discarded.push_back(Set.back());
					Set.pop_back();				n_eliminated++;
				}

			}
			else if (type == CPX_MAX) {
				if (aa == 1) { //it1 = v >  it2
					Discarded.push_back(Set.back());
					Set.pop_back();				n_eliminated++;
				}
				else if (aa == 2) { //it1 = v >_lex it2. They are ordered

				}
				else if (aa == 0) { //The two vectors are equal   it1 = it2
					Discarded.push_back(Set.front());
					Set.pop_front();
					return -1;
				}
				else if (aa == -2) {//it1 = v <_lex it2.
					Set.push_back(*it1);				Set.pop_front();
				}
				else if (aa == -1) { //it1 <  it2. We eliminate it1
					Discarded.push_back(Set.front());
					Set.pop_front();
					return -1;
				}
			}
		}
		else if (Set.size() > 2) { //The list has at least two elements to compare with
			std::advance(it2, 1);
			int aa;

			if (type == CPX_MIN) {
				while (!finish) {
					aa = utils::compare_lex(&(*it1), &(*it2), 0, p);
					if (aa == 1) { //v1 >  v2
						Discarded.push_back(*it1);
						it1 = Set.erase(it1);
						return -1;
					}
					else if (aa == 2) { //If v1 >_lex v2, go to the next element of the list
						++it2;
					}
					else if (aa == 0) {
						Discarded.push_back(*it1);
						it1 = Set.erase(it1);
						return -1;
					}
					else if (aa == -1) { //v1 < v2
						Discarded.push_back(*it2);
						it2 = Set.erase(it2);
						n_eliminated++;
						if (it2 == Set.end()) {
							Set.insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
							Set.pop_front();
							return (n_eliminated);
						}
					}
					else { //Position found
						finish = true;
						bool finish2 = false;

						std::list<point>::iterator it3 = it2;
						while (!finish2) {			//We look for more dominated points
							aa = utils::compare_lex(&(*it1), &(*it3), 0, p);
							if (aa == -1) {
								Discarded.push_back(*it3);
								it3 = Set.erase(it3);
								n_eliminated++;
							}
							else {
								++it3;
							}

							if (it3 == Set.end()) {
								finish2 = true;
							}
						}
					}
					if (it2 == Set.end())
						finish = true;
				}
				Set.insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
				Set.pop_front();
			}
			else if (type == CPX_MAX) {
				while (!finish) {
					aa = utils::compare_lex(&(*it1), &(*it2), 0, p);

					if (aa == 1) { //v1 >  v2
						Discarded.push_back(*it2);
						it2 = Set.erase(it2);
						n_eliminated++;
						if (it2 == Set.end()) {
							Set.insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
							Set.pop_front();
							return (n_eliminated);
						}
					}
					else if (aa == 2) { //Position found
						finish = true;
						bool finish2 = false;

						std::list<point>::iterator it3 = it2;
						while (!finish2) {			//We look for more dominated points
							aa = utils::compare_lex(&(*it1), &(*it3), 0, p);
							if (aa == 1) {
								Discarded.push_back(*it3);
								it3 = Set.erase(it3);
								n_eliminated++;
							}
							else {
								++it3;
							}

							if (it3 == Set.end()) {
								finish2 = true;
							}
						}
					}
					else if (aa == 0) {
						Discarded.push_back(*it1);
						it1 = Set.erase(it1);
						return -1;
					}
					else if (aa == -1) { //v1 < v2
						Discarded.push_back(*it1);
						it1 = Set.erase(it1);
						return -1;
					}
					else { //If v1 >_lex v2, go to the next element of the list
						++it2;
					}
					if (it2 == Set.end())
						finish = true;
				}
				Set.insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
				Set.pop_front();
			}
		}
		else { //PF->size == 1
		}
		return (n_eliminated);
	}

	static void filter_and_order(std::list<std::vector<double>> &Set, std::list<std::vector<double>> &Discarded, const int order, int p) {
		//Create a copy of vectors Set and Discarded (Set2 and Discarded2). Clear Set and Discarded. Insert every element of Set2 into Set avoiding redundancy. New Discarded is created. Add to Discarded the previous elements in Discarded2.
		//order = 1. Set the points with <_lex
		//order = -1. Set the points with >_lex
		int out;
		std::list<std::vector<double>> Set2 = Set;
		std::list<std::vector<double>> Discarded2 = Discarded;
		
		Set.clear();
		Discarded.clear();
		
		for (std::list<std::vector<double>>::iterator it = Set2.begin(); it != Set2.end(); ++it) {
			out = insert_without_domination(order, &(*it), Set, Discarded, p);
		}
		for (std::list<std::vector<double>>::iterator it = Discarded2.begin(); it != Discarded2.end(); ++it) {
			Discarded.push_front(*it);
		}
	}
	
	static bool compare_sets(std::list<std::vector<double>> &A, std::list<std::vector<double>> &B, int p){

		std::list<std::vector<double>> Equal, AminusB, BminusA, DominatedA, DominatedB;
		std::list<std::vector<double>>::iterator it1, it2;
		int value;
		bool equal = true;

		it1 = A.begin();
		it2 = B.begin();

		while ( (it1 != A.end()) && (it2 != B.end()) ){
			value = compare_lex(&(*it1), &(*it2), 0 , p);

			switch(value){
				case -2:
					AminusB.push_back(*it1);
					equal =false;
					++it1;
					break;
				case -1:
					DominatedA.push_back(*it1);
					equal =false;
					++it1;
					break;
				case 0:
					Equal.push_back(*it1);
					++it1; ++it2;
					break;
				case 1:
					DominatedB.push_back(*it2);
					equal =false;
					++it2;
					break;
				case 2:
					BminusA.push_back(*it2);
					equal =false;
					++it2;
					break;
			}		
		}
		while (it1 != A.end()){
			AminusB.push_back(*it1);
			++it1;
		}
		while (it2 != B.end()){
			BminusA.push_back(*it2);
			++it2;
		}

		if (!equal){
			int sksks =3;
		}


		return equal;

	}
	~utils() {}
};

class Box {
public:
#ifdef ON_SCREEN
	int id;						//Box id 
#endif
	point lb;					//Lower bound of the box
	point ub;					//Upper bound of the box
	double value;				//Value of the box
	int direction;				//Position of the box refering to its parent
	Box* parent;				//Parent of the box
	

public:
	Box() {}
	void set(point lb_, point ub_, int direction_, const double *scaling, Box* parent_ = NULL) {
#ifdef ON_SCREEN
		active_boxes++;		
		id = box_id++;
#endif
		lb = lb_;
		ub = ub_;
		value = get_value(lb, ub, scaling);
		direction = direction_;
		parent = parent_;
		
	}
	
private:
	static int compare_lex(Box &v1, Box &v2, int index, int p) {
		
		return (utils::compare_lex(&v1.ub, &v2.ub, 0 , p));
	}
	static double get_volume(const point &lb, const point &ub) {
		//Get the volume of the box
		double volume = 1.0;
		for (int i = 0; i < lb.size(); i++) {
			volume = volume * (ub[i] - lb[i]);
		}
		return volume;
	}
	static double get_scaled_volume(const point &lb, const point &ub, const double *scaling) {
		//Get scaled volume of the box
		double volume = 1.0;
		for (int i = 0; i < lb.size(); i++) {
			volume = volume * ((ub[i] - lb[i]) / scaling[i]);
		}
		return volume;
	}
public:
	int insert_without_domination(std::list<Box *> *L, std::list<Box *> *Eliminated_vectors, const int p) {
		//Insert a new box v in lexicographic order (>_lex) into list L. Dominated vectors are saved in 'Eliminated_vectors'.
		//Return 0 if success
		//Return -1 if its upper bound if dominated by another upper bound (thus, not included in the list)
		//Return a >= 1 , where a is the number of boxes which are deleted from the list, because they are dominated by the new box

		//Outputs of compare_lex.
		//Compare two upper bounds for boxes v1 and v2 in lexicographic order.
		//Return -2 if v1 <_lex v2 and not v1 < v2 ; Return -1 if v1 <v2; Return 0 if v1 = v2; Return 1  if v1 > v2; Return 2 if v1 >_lex v2 and not v1 > v2.  Dominated vectors are saved in 'Eliminated_vectors'

		std::list<Box *>::iterator it1, it2, aux;
		int n_eliminated = 0;
		bool finish = false;

		//Introduce the new point in the front of the list. Then compare it with every element on the list, one by one until its lexicographical position is found.
		L->push_front(this);
		it1 = it2 = L->begin();

		if (L->size() == 2) { //There is only one element to compare with
			std::advance(it2, 1);
			int aa = compare_lex(**it1, **it2, 0, p);

			if (aa == 1) { //it1 = v >  it2
				Eliminated_vectors->push_back(L->back());
				L->pop_back();				n_eliminated++;
			}
			else if (aa == 2) { //it1 = v >_lex it2. They are ordered

			}
			else if (aa == 0) { //The two vectors are equal   it1 = it2
				Eliminated_vectors->push_back(L->front());
				L->pop_front();
				return -1;
			}
			else if (aa == -2) {//it1 = v <_lex it2. 
				L->push_back(*it1);				L->pop_front();
			}
			else if (aa == -1) { //it1 <  it2. We eliminate it1
				Eliminated_vectors->push_back(L->front());
				L->pop_front();
				return -1;
			}

		}
		else if (L->size() > 2) { //The list has at least two elements to compare with
			std::advance(it2, 1);
			int aa;

			while (!finish) {
				aa = compare_lex(**it1, **it2, 0, p);

				if (aa == 1) { //v1 >  v2
					Eliminated_vectors->push_back(*it2);
					it2 = L->erase(it2);
					n_eliminated++;
					if (it2 == L->end()) {
						L->insert(it2, this); //The position of vector v is found. Insert it and remove the front of the list
						L->pop_front();
						return (n_eliminated);
					}
				}
				else if (aa == 2) { //Position found 
					finish = true;
					bool finish2 = false;

					std::list<Box *>::iterator it3 = it2;
					while (!finish2) {			//We look for more dominated points
						aa = compare_lex(**it1, **it3, 0, p);
						if (aa == 1) {
							Eliminated_vectors->push_back(*it3);
							it3 = L->erase(it3);
							n_eliminated++;
						}
						else {
							++it3;
						}
						if (it3 == L->end()) {
							finish2 = true;
						}
					}
				}
				else if (aa == 0) {
					Eliminated_vectors->push_back(*it1);
					it1 = L->erase(it1);
					return -1;
				}
				else if (aa == -1) { //v1 < v2
					Eliminated_vectors->push_back(*it1);
					it1 = L->erase(it1);
					return -1;
				}
				else { //If v1 >_lex v2, go to the next element of the list
					++it2;
				}
				if (it2 == L->end())
					finish = true;
			}
			L->insert(it2, this); //The position of vector v is found. Insert it and remove the front of the list
			L->pop_front();
		}
		else { //PF->size == 1
		}

		return (n_eliminated);
	}
	static void update_to_min_lb(point &LB, const point &lb2) {
		//Point LB is modified.   LB_i = min (LB_i, lb2_i)
		int i;
		
		for (i = 0; i < LB.size(); i++) {
			if (lb2.at(i) < LB.at(i)) {
				LB.at(i) = lb2.at(i);
			}
		}
	}
	int insert_without_domination_and_join(std::list<Box *> *L, std::list<Box *> *Eliminated_vectors, const int p) {
		//Insert a new box v in lexicographic order (>_lex) into the list L, using P-PARTITION. That means that we call to 'join' boxes (update_to_min_lb). 
		//Return 0 if success
		//Return -1 if its upper bound if dominated by another upper bound (thus, not included in the list)
		//Return a >= 1 , where a is the number of boxes which are deleted from the list, because they are dominated by the new box

		//Outputs of compare_lex.
		//Compare two upper bounds for boxes v1 and v2 in lexicographic order.
		//Return -2 if v1 <_lex v2 and not v1 < v2 ; Return -1 if v1 <v2; Return 0 if v1 = v2; Return 1  if v1 > v2; Return 2 if v1 >_lex v2 and not v1 > v2.  Dominated vectors are saved in 'Eliminated_vectors'
		std::list<Box *>::iterator it1, it2, aux;
		int n_eliminated = 0;
		bool finish = false;

		//Introduce the new point in the front of the list. Then compare it with every element on the list, one by one until its lexicographical position is found.
		L->push_front(this);
		it1 = it2 = L->begin();

		
		if (L->size() == 2) { //There is only one element to compare with
			std::advance(it2, 1);
			int aa = compare_lex(**it1, **it2, 0, p);

			if (aa == 1) { //it1 = v >  it2
				update_to_min_lb((*it1)->lb, (*it2)->lb);
				Eliminated_vectors->push_back(L->back());
				L->pop_back();				n_eliminated++;
			}
			else if (aa == 2) { //it1 = v >_lex it2. They are ordered

			}
			else if (aa == 0) { //The two vectors are equal   it1 = it2
				update_to_min_lb((*it2)->lb, (*it1)->lb);
				Eliminated_vectors->push_back(L->front());		////
				L->pop_front();
				return -1;
			}
			else if (aa == -2) {//it1 = v <_lex it2. 
				L->push_back(*it1);				L->pop_front();
			}
			else if (aa == -1) { //it1 <  it2. We eliminate it1
				update_to_min_lb((*it2)->lb, (*it1)->lb);
				Eliminated_vectors->push_back(L->front());		////
				L->pop_front();
				return -1;
			}

		}
		else if (L->size() > 2) { //The list has at least two elements to compare with
			std::advance(it2, 1);
			int aa;

			while (!finish) {
				aa = compare_lex(**it1, **it2, 0, p);

				if (aa == 1) { //v1 >  v2
					update_to_min_lb((*it1)->lb, (*it2)->lb);
					Eliminated_vectors->push_back(*it2);
					it2 = L->erase(it2);
					n_eliminated++;
					if (it2 == L->end()) {
						L->insert(it2, this); //The position of vector v is found. Insert it and remove the front of the list
						L->pop_front();
						return (n_eliminated);
					}
				}
				else if (aa == 2) { //Position found 
					finish = true;
					bool finish2 = false;

					std::list<Box *>::iterator it3 = it2;
					while (!finish2) {			//We look for more dominated points
						aa = compare_lex(**it1, **it3, 0, p);
						if (aa == 1) {
							update_to_min_lb((*it1)->lb, (*it3)->lb);
							Eliminated_vectors->push_back(*it3);
							it3 = L->erase(it3);
							n_eliminated++;
						}
						else {
							++it3;
						}

						if (it3 == L->end()) {
							finish2 = true;
						}
					}
				}
				else if (aa == 0) {
					update_to_min_lb((*it2)->lb, (*it1)->lb);
					Eliminated_vectors->push_back(*it1);		////
					it1 = L->erase(it1);
					return -1;
				}
				else if (aa == -1) { //v1 < v2
					update_to_min_lb((*it2)->lb, (*it1)->lb);
					Eliminated_vectors->push_back(*it1);		////
					it1 = L->erase(it1);
					return -1;
				}
				else { //If v1 >_lex v2, go to the next element of the list
					++it2;
				}
				if (it2 == L->end())
					finish = true;
			}
			L->insert(it2, this); //The position of vector v is found. Insert it and remove the front of the list
			L->pop_front();
		}
		else { //PF->size == 1
		}
		
		return (n_eliminated);

	}
	static double get_value(const point &lb, const point &ub, const double *scaling) {
		//Get the value of the box
		if (scaling == NULL) return(Box::get_volume(lb, ub));
		else return(Box::get_scaled_volume(lb, ub, scaling));
	}
	
	~Box() {
#ifdef ON_SCREEN
		active_boxes--;
		printf("\nDestroy box %d", id);
#endif
		
	}
};

class Metrics {
public:
	std::list<std::vector<double>> set;
	std::list<std::vector<double>> PF;
	int ONVG;
	double ONVGR, GD, IGD, HV;

public:
	Metrics() {}

	Metrics(const std::list<std::vector<double>> &thisset, const std::list<std::vector<double>> &thisPF) {
		set = thisset;	PF = thisPF; ONVG = 0; ONVGR = -1; GD = -1; IGD = -1; HV = -1;
	}
private:
	//Euclidean distance
	double distance(const std::vector<double> &a, const std::vector<double> &b) {
		double distance = 0.0;
		for (int i = 0; i < a.size(); i++) {
			distance += pow(a[i] - b[i], 2.0);
		}
		return sqrt(distance);
	} // distance
	  // Gets the distance between a point and the nearest one in a given front 
	double distanceToClosedPoint(const std::vector<double> &point, std::list<std::vector<double>> &front) {
		double minDistance = MAX_DOUBLE;
		double aux;
		for (std::list<std::vector<double>>::iterator it = front.begin(); it != front.end(); ++it) {
			if (minDistance != 0) {
				aux = distance(point, *it);
				if (aux < minDistance) minDistance = aux;
			}
		}
		return minDistance;

	} // distanceToClosedPoint
	  // Calculate GD metric
	double GenerationalDistance(std::list<std::vector<double>> &S, std::list<std::vector<double>> &P) {
		double sum = 0.0;
		double GD = -1;
		if ((S.size() == 0) || (P.size() == 0)) return GD;
		for (std::list<std::vector<double>>::iterator it = S.begin(); it != S.end(); ++it) {
			sum += pow(distanceToClosedPoint(*it, P), 2);
		}
		sum = sqrt(sum);
		GD = sum / S.size();
		return GD;
	}
	//Calculate IGD metric
	double InvertedGenerationalDistance(std::list<std::vector<double>> &S, std::list<std::vector<double>> &P) {
		return (GenerationalDistance(P, S));
	}
#ifdef HV_WFG
private:
	double Calculate_Hypervolume_WFG2(std::list<std::vector<double>> *PF, std::vector<double> *reference_point, int p, int type) {
		FILE *fp;
		int i;
		double hypervolume = 0;

		std::string file = "_aux_";

		//Convert referente point to a pointer of a vector
		double *ref_point = (double *)malloc(p * sizeof(double));
		for (i = 0; i < p; i++) ref_point[i] = reference_point->at(i);

		fp = fopen(file.c_str(), "w");
		fprintf(fp, "#");

		if (type == CPX_MAX) { //In MAX problems we create a new file of type MIN
			for (std::list<std::vector<double>>::iterator it = PF->begin(); it != PF->end(); ++it) {
				fprintf(fp, "\n");
				for (i = 0; i < p; i++) {
					fprintf(fp, "%lf ", -(*it).at(i));
				}
			}
			for (i = 0; i < p; i++) ref_point[i] = -ref_point[i];
		}
		else {
			for (std::list<std::vector<double>>::iterator it = PF->begin(); it != PF->end(); ++it) {
				fprintf(fp, "\n");
				for (i = 0; i < p; i++) {
					fprintf(fp, "%lf ", (*it).at(i));
				}
			}
		}
		fprintf(fp, "\n#");
		fclose(fp);

		//Call to While & Bradstree & Barone software (WFG)
		char *archive2 = strcpy((char*)malloc(file.length() + 1), file.c_str());
		hypervolume = WFG_HV(archive2, ref_point);

		remove(archive2);
		free(archive2);
		free(ref_point);

		return (hypervolume);
	}
public:
	void get_HV(int delta, int &type, int &p) {
		std::vector<double> ref_point(p);
		get_reference_point(ref_point, delta, type);
		HV = Calculate_Hypervolume_WFG2(&set, &ref_point, p, type);
	}
	void get_HV(std::vector<double> &ref_point, int &type, int &p) {
		HV = Calculate_Hypervolume_WFG2(&set, &ref_point, p, type);
	}
#endif

#ifdef HV_FONSECA
	private:
		double Calculate_Hypervolume_Fonseca(std::list<std::vector<double>> *PF, std::vector<double> *reference_point, int type, int p) {

		
		FILE *fp;
		int i;
		double hypervolume = 0;

		std::string file = "_aux_";

		//Convert referente point to a pointer of a vector
		double *ref_point = (double *)malloc(p * sizeof(double));
		for (i = 0; i < p; i++) ref_point[i] = reference_point->at(i);

		fp = fopen(file.c_str(), "w");
		

		if (type == CPX_MAX) { //In MAX problems we create a new file of type MIN
			for (std::list<std::vector<double>>::iterator it = PF->begin(); it != PF->end(); ++it) {
				fprintf(fp, "\n");
				for (i = 0; i < p; i++) {
					fprintf(fp, "%lf ", -(*it).at(i));
				}
			}
			for (i = 0; i < p; i++) ref_point[i] = -ref_point[i];
		}
		else {
			for (std::list<std::vector<double>>::iterator it = PF->begin(); it != PF->end(); ++it) {
				fprintf(fp, "\n");
				for (i = 0; i < p; i++) {
					fprintf(fp, "%lf ", (*it).at(i));
				}
			}
		}
		fprintf(fp, "\n");
		fclose(fp);



		
		//Call to Fonseca & Lopez-Ibanez & Paquete software 
		char *archive2 = strcpy((char*)malloc(file.length() + 1), file.c_str());
		//hypervolume = hv_file(archive2, ref_point, NULL, NULL, &p); //FONSECA - MANUEL LOPEZ-IBA�EZ


		remove(archive2);
		if (archive2 != NULL) free(archive2);
		if (ref_point != NULL) free(ref_point);

		return (hypervolume);
	}
	public:
		void get_HVf(int delta, int &type, int &p) {
		std::vector<double> ref_point(p);
		get_reference_point(ref_point, delta, type);
		HV = Calculate_Hypervolume_Fonseca(&set, &ref_point, type, p);
	}
	void get_HVf(std::vector<double> &ref_point, int &type, int &p) {
		HV = Calculate_Hypervolume_Fonseca(&set, &ref_point, type, p);
	}
#endif

	void get_reference_point(std::vector<double> &ref_point, int delta, int type) {
		int p = (int)ref_point.size();
		std::vector<double> Nadir(p);
		if (type == CPX_MAX) {
			//Calculate Nadir point and returns Nadir - delta
			for (int i = 0; i < p; i++) Nadir.at(i) = INT_MAX;
			for (std::list<std::vector<double>>::iterator it = PF.begin(); it != PF.end(); ++it) {
				for (int i = 0; i < p; i++) {
					if ((*it).at(i) < Nadir.at(i)) Nadir.at(i) = (*it).at(i);
				}
			}
			for (int i = 0; i < p; i++) {
				ref_point.at(i) = Nadir.at(i) - delta;
			}
		}
		else if (type == CPX_MIN) {
			//Calculate Nadir point and returns Nadir + delta
			for (int i = 0; i < p; i++) Nadir.at(i) = -INT_MAX;
			for (std::list<std::vector<double>>::iterator it = PF.begin(); it != PF.end(); ++it) {
				for (int i = 0; i < p; i++) {
					if ((*it).at(i) > Nadir.at(i)) Nadir.at(i) = (*it).at(i);
				}
			}
			for (int i = 0; i < p; i++) {
				ref_point.at(i) = Nadir.at(i) + delta;
			}
		}
	}

public:
	bool get_set(std::list<solution> &SET) {
		std::list<solution>::iterator it = SET.begin();
		int p = int((*it).z.size());
		while (it != SET.end()) {
			std::vector<double> newvector(p);
			for (int j = 0; j < p; j++) newvector[j] = (*it).z[j];
			set.push_back(newvector);
			++it;
		}
		return true;
	}
	bool get_pf(std::string location, int p) {
		int i;

		double x;
		int count;

		FILE *fp2;

		fp2 = fopen(location.c_str(), "r");
		if (fp2 == NULL) {
			printf("\nArchive %s not found.", location.c_str());
			return false;
		}
		else {
			count = 0;
			std::vector<double> newv(p);
			while (!feof(fp2)) {
				for (i = 0; i < p; i++) {
					fscanf(fp2, "%lf", &x);
					newv.at(i) = x;
				}
				PF.push_back(newv);
				count++;
			}
			fclose(fp2);
		}
		return true;
	}
	void get_GD() {
		GD = GenerationalDistance(set, PF);
	}
	void get_IGD() {
		IGD = GenerationalDistance(PF, set);
	}
	void get_ONVG() {
		ONVG = (int(set.size()));
	}
	void get_ONVGR() {
		if (PF.size() == 0) ONVGR = -1;
		ONVGR = double(set.size()) / double(PF.size());
	}


public:
	~Metrics() {
		//set.clear(); PF.clear();
	} //Destructor
};

class Compare_value {
public:
	bool operator()(Box *x, Box *y) const {
		return x->value > y->value;
	}
};

class CPLEX {
public:
	CPXENVptr *env;
	CPXLPptr *lp;
	int n_iterations = 0;
	double t_solver = 0;
	int *indexes;
	double t_max;
	int max_iter;
	double tilim;
private:
	int beg = 0, effort = 0; 
public:
	CPLEX() {
		env = new(CPXENVptr); //Create CPLEX environment
		lp = new(CPXLPptr);	//Create associated problem to CPLEX environment
	}
	void create(std::string &name) {
		int status = 0;
		const char *namec = name.c_str();
		*env = CPXopenCPLEX(&status);
		*lp = CPXcreateprob(*env, &status, namec);
	}
	void set_parameters(const int PARALLEL_, const int THREADS_, double EPAGAP_ = 1e-6, double EPINT_ = 1e-5) {
		int status = 0;
		status = CPXsetintparam(*env, CPX_PARAM_PARALLELMODE, PARALLEL); //-1 CPX_PARALLEL_OPPORTUNISTIC;   0 (default) CPX_PARALLEL_AUTO  1 CPX_PARALLEL_DETERMINISTIC
		status = CPXsetintparam(*env, CPX_PARAM_THREADS, THREADS);	//0 (default, CPLEX decides) ;	1 (sequential, single Thread) ;		N (uses up to N threads)
		status = CPXsetdblparam(*env, CPX_PARAM_EPAGAP, EPAGAP); //Default 1e-06. Sets an absolute tolerance on the gap between the best integer objective and the objective of the best node remaining.
		status = CPXsetdblparam(*env, CPX_PARAM_EPINT, EPINT); //Default 1e-05. CPLEX tolerance for integer variables
	}
	void set_other(int nvar, double &t_max_, int &max_iter_ ) {
		indexes = (int *)malloc(nvar * sizeof(int));
		for (int i = 0; i < nvar; i++) indexes[i] = i;
		t_max = t_max_;
		max_iter = max_iter_;
	}
	double LP_Solve_z(int &stat, CPXLPptr *lp_) {
		//Solve a LP problem, returning the stat and objective value
		TIEMPO t;
		t.init();
		double z;
		CPXlpopt(*env, *lp_);
		CPXgetobjval(*env, *lp_, &z);
		stat = CPXgetstat(*env, *lp_);
		n_iterations++;
		t.acum();
		t_solver += t.value();
		return z;
	}
	double ILP_Solve_z(int &stat, double **x, int &nvar) {
		TIEMPO t;
		t.init();
		double z;
		//Solve a MIP problem, returning the stat and objective value
		CPXmipopt(*env, *lp);
		CPXgetobjval(*env, *lp, &z);
		stat = CPXgetstat(*env, *lp);
		CPXgetx(*env, *lp, *x, 0, nvar - 1);
		t.acum();
		t_solver += t.value();
		n_iterations++;
		return z;
	}
	void set_time_for_solver(TIEMPO &t1) {
		//Set the maximum time available for CPLEX, which is t_max - t1. Return the previous time limit for CPLEX
		t1.acum();
		double currenttime = t_max - t1.value();
		if (currenttime < tilim) {
			currenttime = std::max(1e-74, currenttime);
			CPXsetdblparam(*env, CPX_PARAM_TILIM, currenttime);
		}
	}
	/*
	void assign_mip_start(Box &B, int &nvar) {
		int nummipstart = CPXgetnummipstarts(*env, *lp);
		if (nummipstart > 0) CPXdelmipstarts(*env, *lp, 0, nummipstart - 1);
		CPXaddmipstarts(*env, *lp, 1, nvar, &beg, indexes, B.last_solution, &effort, NULL);
		//string abc = "hh" + to_string(B->box_id) + ".lp";
		//CPXwritemipstarts(*P->env, *P->lp, abc.c_str(), 0, CPXgetnummipstarts(*P->env, *P->lp) - 1);
	}
	*/
	~CPLEX() {
		free(indexes);
		CPXcloseCPLEX(env);
	}
};

class Solution {
public:
	std::list<solution> Set;
	std::list<solution> Discarded;
	point Bound_I, Bound_N;
	double *scaling;
	int dim;
public:
	Solution() {}

private:
	int compare_lex(solution *v1, solution *v2, int index, int &p) {

		return (utils::compare_lex(&v1->z, &v2->z , index, p));
	
	}
public:
	void copy_set_from_file(std::string fichero, int p) {
	//Clear Set and Discarded, and read new Set from file
		Set.clear();
		Discarded.clear();
		
		int i;
		FILE *fp;
		fp = fopen(fichero.c_str(), "r");

		while (!feof(fp)) {
			std::vector<double> v(p);
			solution NS;
			for (i = 0; i < p; i++) {
				fscanf(fp, "%lf", &v[i]);
			}
			NS.time_point = 0;
			NS.z = v;
			Set.push_back(NS);
		}
		//printf("\nCargado PF de nombre %s con %d elementos", fichero.c_str(), int(L->size()));

	}
	int insert_without_domination(int type, solution *v, int &p) {
		//Insert a new nondominated point z in lexicographic order into the list L.
		//If type = 1, use <_lex		//If type = -1, use >_lex
		//Return 0 if success
		//Return -1 if the solution is dominated (thus, not included in the list)
		//Return a >= 1 , where a is the number of solutions which are deleted from the list, because they are dominated by the new solution

		//Outputs of compare_lex. //Compare two p-dimensional vectors v1 and v2 in lexicographic order.
		//Return -2 if v1 < v2; Return -1 if v1 <_lex v2 and not v1 < v2; Return 0 if v1 = v2; Return 1  if v1 > v2; Return 2 if v1 >_lex v2 and not v1 > v2

		std::list<solution>::iterator it1, it2, aux;
		int n_eliminated = 0;
		bool finish = false;

		//Introduce the new point in the front of the list. Then compare it with every element on the list, one by one until its lexicographical position is found.
		Set.push_front(*v);
		it1 = it2 = Set.begin();

		if (Set.size() == 2) { //There is only one element to compare with
			std::advance(it2, 1);
			int aa = Solution::compare_lex(&(*it1), &(*it2), 0, p);

			if (type == CPX_MIN) {
				if (aa == 1) { //it1 = v >  it2
					Discarded.push_back(Set.front());
					Set.pop_front();
					return -1;
				}
				else if (aa == 2) { //it1 = v >_lex it2
					Set.push_back(*it1);				Set.pop_front();
				}
				else if (aa == 0) { //The two vectors are equal   it1 = it2
					Discarded.push_back(Set.front());
					Set.pop_front();
					return -1;
				}
				else if (aa == -2) {//it1 = v <_lex it2. They are ordered

				}
				else if (aa == -1) { //it1 <  it2. We eliminate it2
					Discarded.push_back(Set.back());
					Set.pop_back();				n_eliminated++;
				}

			}
			else if (type == CPX_MAX) {
				if (aa == 1) { //it1 = v >  it2
					Discarded.push_back(Set.back());
					Set.pop_back();				n_eliminated++;
				}
				else if (aa == 2) { //it1 = v >_lex it2. They are ordered

				}
				else if (aa == 0) { //The two vectors are equal   it1 = it2
					Discarded.push_back(Set.front());
					Set.pop_front();
					return -1;
				}
				else if (aa == -2) {//it1 = v <_lex it2.
					Set.push_back(*it1);				Set.pop_front();
				}
				else if (aa == -1) { //it1 <  it2. We eliminate it1
					Discarded.push_back(Set.front());
					Set.pop_front();
					return -1;
				}
			}

		}
		else if (Set.size() > 2) { //The list has at least two elements to compare with
			std::advance(it2, 1);
			int aa;


			if (type == CPX_MIN) {
				while (!finish) {
					aa = Solution::compare_lex(&(*it1), &(*it2), 0, p);
					if (aa == 1) { //v1 >  v2
						Discarded.push_back(*it1);
						it1 = Set.erase(it1);
						return -1;
					}
					else if (aa == 2) { //If v1 >_lex v2, go to the next element of the list
						++it2;
					}
					else if (aa == 0) {
						Discarded.push_back(*it1);
						it1 = Set.erase(it1);
						return -1;
					}
					else if (aa == -1) { //v1 < v2
						Discarded.push_back(*it2);
						it2 = Set.erase(it2);
						n_eliminated++;
						if (it2 == Set.end()) {
							Set.insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
							Set.pop_front();
							return (n_eliminated);
						}
					}
					else { //Position found
						finish = true;
						bool finish2 = false;

						std::list<solution>::iterator it3 = it2;
						while (!finish2) {			//We look for more dominated points
							aa = Solution::compare_lex(&(*it1), &(*it3), 0, p);
							if (aa == -1) {
								Discarded.push_back(*it3);
								it3 = Set.erase(it3);
								n_eliminated++;
							}
							else {
								++it3;
							}

							if (it3 == Set.end()) {
								finish2 = true;
							}
						}
					}
					if (it2 == Set.end())
						finish = true;
				}

				Set.insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
				Set.pop_front();
			}
			else if (type == CPX_MAX) {
				while (!finish) {
					aa = Solution::compare_lex(&(*it1), &(*it2), 0, p);

					if (aa == 1) { //v1 >  v2
						Discarded.push_back(*it2);
						it2 = Set.erase(it2);
						n_eliminated++;
						if (it2 == Set.end()) {
							Set.insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
							Set.pop_front();
							return (n_eliminated);
						}
					}
					else if (aa == 2) { //Position found
						finish = true;
						bool finish2 = false;

						std::list<solution>::iterator it3 = it2;
						while (!finish2) {			//We look for more dominated points
							aa = Solution::compare_lex(&(*it1), &(*it3), 0, p);
							if (aa == 1) {
								Discarded.push_back(*it3);
								it3 = Set.erase(it3);
								n_eliminated++;
							}
							else {
								++it3;
							}

							if (it3 == Set.end()) {
								finish2 = true;
							}
						}
					}
					else if (aa == 0) {
						Discarded.push_back(*it1);
						it1 = Set.erase(it1);
						return -1;
					}
					else if (aa == -1) { //v1 < v2
						Discarded.push_back(*it1);
						it1 = Set.erase(it1);
						return -1;
					}
					else { //If v1 >_lex v2, go to the next element of the list
						++it2;
					}
					if (it2 == Set.end())
						finish = true;
				}
				Set.insert(it2, *v); //The position of vector v is found. Insert it and remove the front of the list
				Set.pop_front();
			}
		}
		else { //PF->size == 1
		}

		return (n_eliminated);

	}
	bool check(int &p) {
		//Check if the Set have non-dominated points
		Solution Set2;
		int out;
		for (std::list<solution>::iterator it = Set.begin(); it != Set.end(); ++it) {
			out = Set2.insert_without_domination(1, &(*it), p);
			if (out != 0) {
				return false;
			}
		}
		return true;
	}
	void filter(const int order, int p) {
		//Create a copy of Set and Discarded (Set2 and Discarded2). Clear Set and Discarded. Insert every element of Set2 into Set avoiding redundancy. New Discarded is created. Add to Discarded the previous elements in Discarded2.
		//order = 1. Set the points with <_lex
		//order = -1. Set the points with >_lex
		int out;
		std::list<solution> Set2 = Set;
		std::list<solution> Discarded2 = Discarded;
		
		Set.clear();
		Discarded.clear();
		
		for (std::list<solution>::iterator it = Set2.begin(); it != Set2.end(); ++it) {
			out = insert_without_domination(order, &(*it),p);
		}
		for (std::list<solution>::iterator it = Discarded2.begin(); it != Discarded2.end(); ++it) {
			Discarded.push_front(*it);
		}
	}
	void copy_set_to_file(std::string file_name, int lenght) {
		FILE *fp;
		int i, j;

		fp = fopen((file_name).c_str(), "w");

		if (Set.size() > 0) {
			std::list<solution>::iterator it = Set.begin();
			for (j = 0; j < Set.size() - 1; j++) {
				for (i = 0; i < lenght - 1; i++) fprintf(fp, "%1.0f ", it->z[i]);
				fprintf(fp, "%1.0f\n", it->z[i]);
				++it;
			}
			for (i = 0; i < lenght - 1; i++) fprintf(fp, "%1.0f ", it->z[i]);
			fprintf(fp, "%1.0f", it->z[i]);
		}
		fclose(fp);
	}
	void copy_set_and_time_to_file(std::string file_name, int lenght) {
		FILE *fp;
		int i, j;

		fp = fopen((file_name).c_str(), "w");

		if (Set.size() > 0) {
			std::list<solution>::iterator it = Set.begin();
			for (j = 0; j < Set.size() - 1; j++) {
				for (i = 0; i < lenght - 1; i++) fprintf(fp, "%.0f ", it->z[i]);
				fprintf(fp, "%.0f %.4f\n", it->z[i], it->time_point);
				++it;
			}
			for (i = 0; i < lenght - 1; i++) fprintf(fp, "%1.0f ", it->z[i]);
			fprintf(fp, "%.0f %.4f", it->z[i], it->time_point);
		}
		fclose(fp);
	}
	~Solution() {
	}
};

class MOILP { //Multi-objective integer linear problem
public:
	CPLEX cplex;						//CPLEX object
	int dimension;						//Dimension
	int n_var;							//Number of decision variables
	int n_const;						//Number of constraints
	double **F;							//Matrix of objective coefficients
	Solution APF;						//Solution sets
	Stat stat;							//Solution stats
	int original_sense;					//1 MIN;  -1 MAX
	int counter;

public:
	MOILP() {}
	void set_model(std::string &path, std::string &name, double tmax, int maxcplexiter, int PARALLEL_, int THREADS_, double TILIM_) {

		counter = 0;

		cplex.create(name);
		cplex.set_parameters(PARALLEL_, THREADS_, TILIM_);
		cplex.tilim = TILIM_;

		//Read .lp CPLEX model and objective archive. If the problem is MAX, convert it into MIN
		FILE *fp;
		int i, j, status;
		double bound;
		std::string str = path + ".lp";

		//Read CPLEX problem
		status = CPXreadcopyprob(*cplex.env, *cplex.lp, str.c_str(), NULL);
		if (status != 0) utils::error_CPLEX_status(status);

		str = path + ".obj";

		fp = fopen(str.c_str(), "r");
		if (fp == NULL) utils::error_open_file(str.c_str());

		fscanf(fp, "%d", &dimension); //Read number of objectives
		fscanf(fp, "%d", &n_var); //Read number of variables
		fscanf(fp, "%d", &n_const); //Read number of constraints

		APF.scaling = (double *)malloc(dimension * sizeof(double));
		

		//Creating the objective costs matrix
		F = (double **)malloc(dimension * sizeof(double *));
		for (i = 0; i < dimension; i++) {
			F[i] = (double *)malloc(n_var * sizeof(double));
			for (j = 0; j < n_var; j++) {
				fscanf(fp, "%lf", &F[i][j]);	//Read costs for every objective
			}
		}
		fclose(fp);

		//MAX problems are transformed into MIN problems
		original_sense = CPXgetobjsen(*cplex.env, *cplex.lp);
		if (original_sense == CPX_MAX) {
			for (i = 0; i < dimension; i++) {
				for (j = 0; j < n_var; j++) {
					F[i][j] = -F[i][j];
				}
			}
			CPXchgobjsen(*cplex.env, *cplex.lp, CPX_MIN);
		}

		//Max time and point limit
		if (tmax == 0) tmax =MAX_DOUBLE;
		if (maxcplexiter == 0) maxcplexiter = MAX_INT;

		cplex.set_other(n_var, tmax, maxcplexiter);

		//Bound for Nadir and Ideal point	
		bound = (original_sense == CPX_MIN) ? -CPX_INFBOUND : CPX_INFBOUND;
		APF.Bound_I.resize(dimension);		APF.Bound_N.resize(dimension);
		for (i = 0; i < dimension; i++) {
			APF.Bound_I[i] = bound;
			APF.Bound_N[i] = -bound;
		}
	}
	void create_chalmet_model(void) {
		//Given the original model, create Chalmet et al model
		// MIN (f1 + .... + fp)
		// s.t.		x in X
		//			fi(x) <= ki		for i = 1,...,p
		//
		int i, j, matbeg = 0;
		double *new_c = (double *)malloc(n_var * sizeof(double));

		//The sense of the new constraints are <= because we have a MIN problem
		const char *sense = "L";

		for (i = 0; i < n_var; i++) {
			new_c[i] = 0;
			for (j = 0; j < dimension; j++) {
				new_c[i] += F[j][i];
			}
		}
		CPXchgobj(*cplex.env, *cplex.lp, n_var, cplex.indexes, new_c); //min sum(fi)
		for (i = 0; i < dimension; i++) {
			CPXaddrows(*cplex.env, *cplex.lp, 0, 1, n_var, 0, sense, &matbeg, cplex.indexes, F[i], NULL, NULL); //new constraint
		}
		free(new_c);
	}
	void create_tchebycheff_model(void) {
		//Given the original model, create the Tchebycheff model
		// MIN ((f1 + .... + fp) - alpha)
		// s.t.		x in X
		//			fi(x) <= ki - alpha		for i = 1,...,p
		//			alpha >= 0
		//

		int i, j, matbeg = 0;
		double *new_c = (double *)malloc((n_var + 1) * sizeof(double));

		//The sense of the new constraints are <=, since we have a MIN problem
		const char *sense = "L";

		for (i = 0; i < n_var; i++) {
			new_c[i] = 0;
			for (j = 0; j < dimension; j++) {
				new_c[i] += F[j][i];
			}
		}

		int index_alpha = n_var;
		double coef_alpha = -1.0;
		CPXaddcols(*cplex.env, *cplex.lp, 1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
		CPXchgobj(*cplex.env, *cplex.lp, n_var, cplex.indexes, new_c); //min (eps*sum(fi))
		CPXchgobj(*cplex.env, *cplex.lp, 1, &index_alpha, &coef_alpha); //min (eps*sum(fi) - alpha)

		int numrows = CPXgetnumrows(*cplex.env, *cplex.lp);

		for (i = 0; i < dimension; i++) {
			std::string rowname = "c" + std::to_string(n_const + i + 1);
			char *rownamechar = new char[rowname.length() + 1];
			strcpy(rownamechar, rowname.c_str());
			CPXaddrows(*cplex.env, *cplex.lp, 0, 1, n_var, 0, sense, &matbeg, cplex.indexes, F[i], NULL, &rownamechar); //new constraint   
			numrows++;
			CPXchgcoef(*cplex.env, *cplex.lp, numrows - 1, n_var, 1.0);	//fi(x)+alpha <= 0
			free(rownamechar);
		}
		free(new_c);
	}
	void set_constraint_extreme_box(point &UB) {
		//Set the rhs values for the new constraints (fi(x) <= ki) or (fi(x) - alpha <= ki)  in the current iteration
		//ki = Ni - DELTA, being Ni the bound of the local nadir point (MIN problems)
		//point *N = UB;

		int *indices = (int *)malloc(dimension * sizeof(int));
		double *values = (double *)malloc(dimension * sizeof(double));

		for (int i = 0; i < dimension; i++) {
			indices[i] = n_const + i;
			values[i] = UB[i] - DELTA;	//Subtract a small positive value 
		}
		CPXchgrhs(*cplex.env, *cplex.lp, dimension, indices, values); //s.t.  fi <= ki,  i = 1,...,p 
		free(indices); free(values);
	}
	
	void calculate_image_of_x(double *x, point &z) {
		//Given a vector-solution x, return the point z = (f1(x), ... , fp(x))
		int i, j;
		
		for (i = 0; i < dimension; i++) {
			z[i] = 0.0;
			for (j = 0; j < n_var; j++) {
				z[i] += F[i][j] * x[j];
			}
			if (z[i] - round(z[i]) > DELTA) {
				printf("\nPOINT z EXCEEDS DELTA ERROR");			warning = true;
			}
			z[i] = round(z[i]);
		}
	}
	void transform_pareto_front_if_neccesary(void) {
		//We have used a MIN problem in CPLEX. If the original problem was MAX, we transform the solution
		int i, j;
		if (original_sense == CPX_MAX) {
			for (std::list<solution>::iterator it = APF.Set.begin(); it != APF.Set.end(); ++it) {
				for (j = 0; j < dimension; j++) {
					it->z[j] = -(it->z[j]);
				}
			}
			for (std::list<solution>::iterator it = APF.Discarded.begin(); it != APF.Discarded.end(); ++it) {
				for (j = 0; j < dimension; j++) {
					it->z[j] = -(it->z[j]);
				}
			}
			//Return to the original cost values
			for (i = 0; i < dimension; i++) {
				for (j = 0; j < n_var; j++) {
					F[i][j] = -(F[i][j]);
				}
				APF.Bound_N[i] = -APF.Bound_N[i];
				APF.Bound_I[i] = -APF.Bound_I[i];
			}
		}
	}
	void copy_stats_to_file(std::string file_name, Input &input) {
		FILE *fp;
		fp = fopen(file_name.c_str(), "a+");
		fprintf(fp, "Name: %s", input.name.c_str());
		std::string algorithm = input.partition + "_" + input.model + "_" + input.value + "_" + input.set_of_boxes;
		fprintf(fp, "\nAlgorithm: %s", algorithm.c_str());
		fprintf(fp, "\nTILIM: %.2f", input.TILIM);
		double maxt = input.tmax;
		if (maxt == MAX_DOUBLE) fprintf(fp, "\nMax_time: ---");
		else fprintf(fp, "\nMax_time: %.2f", input.tmax);
		int maxtotalexec = input.maxcplexiter;
		if (maxtotalexec == MAX_INT)	fprintf(fp, "\nMax_total_executions: ---");
		else fprintf(fp, "\nMax_total_executions: %d", input.maxcplexiter);
		fprintf(fp, "\n\nBoundIdeal ( ");
		for (int i = 0; i < dimension; i++) fprintf(fp, "%.1f  ", APF.Bound_I[i]);
		fprintf(fp, ")\tBoundNadir ( ");
		for (int i = 0; i < dimension; i++) fprintf(fp, "%.1f  ", APF.Bound_N[i]);
		fprintf(fp, ")");
		fprintf(fp, "\n|S| = %d", int(APF.Set.size()));
		fprintf(fp, "\nTotal iterations: %d", cplex.n_iterations);
		fprintf(fp, "\nEmpty boxes: %d", stat.empty_boxes);
		fprintf(fp, "\nEliminated boxes because of domination: %d", stat.dominated_boxes);
		fprintf(fp, "\nEliminated solutions because of domination: %d", stat.dominated_points);
		fprintf(fp, "\nRemaining boxes to explore: %d", stat.remaining_boxes);
		double percentage;
		if (stat.total_time == 0) percentage = 0;
		else percentage = round(100 * (cplex.t_solver / stat.total_time));
		fprintf(fp, "\nTime used by the solver: %.2f  (%.0f ./.)", cplex.t_solver, percentage);
		fprintf(fp, "\nTotal time: %.2f (%.0f ./.)", stat.total_time, 100 - percentage);
		fclose(fp);
	}

	~MOILP() {
		for (int i = 0; i < dimension; i++)	free(F[i]);
		free(F);
		if (APF.scaling != NULL) 
			free(APF.scaling);
	}
	public:
		bool estimate_bounds_by_linear_relaxation(TIEMPO &t, std::string &value) {
			//Estimate Ideal point and Nadir point of the problem, by solving the linear relaxation. We always consider MIN problems.
			//In vector scaling we save the range (bound_of_nadir[i] - ideal[i]) for every component
			int i, stat;
			double obj;

			CPXLPptr lp2 = CPXcloneprob(*cplex.env, *cplex.lp, &stat);	//Create a copy of the problem

			stat = CPXchgprobtype(*cplex.env, lp2, CPXPROB_LP);			//Change problem to LP

			int cur_objsense = CPXgetobjsen(*cplex.env, lp2); //CPX_MIN 1	//CPX_MAX -1

			printf("\nEstimating bounds by Linear Relaxation...");

			
			//Ideal point
			for (i = 0; i < dimension; i++) {
				cplex.set_time_for_solver(t);
				stat = CPXchgobj(*cplex.env, lp2, n_var, cplex.indexes, F[i]);		//min(f_i)
				obj = cplex.LP_Solve_z(stat, &lp2);
				//CPXgetx(*cplex.env, lp2, x, 0 , 6);
				APF.Bound_I[i] = floor(obj);
				if (stat != CPX_STAT_OPTIMAL) {
					printf("\nNON OPTIMAL SOLUTION IS FOUND");
					return false;
				}
			}

			//Upper bound for Nadir point
			CPXchgobjsen(*cplex.env, lp2, -cur_objsense);
			for (i = 0; i < dimension; i++) {
				cplex.set_time_for_solver(t);
				stat = CPXchgobj(*cplex.env, lp2, n_var, cplex.indexes, F[i]);		//max(f_i)
				obj = cplex.LP_Solve_z(stat, &lp2);
				//CPXgetx(*cplex.env, lp2, x, 0 , 6);
				APF.Bound_N[i] = ceil(obj + DELTA);
				if (stat != CPX_STAT_OPTIMAL) {
					printf("\nNON OPTIMAL SOLUTION IS FOUND");
					return false;
				}
			}

			//Restore
			CPXchgobjsen(*cplex.env, lp2, cur_objsense);

			//Calculating ranges for every component
			if (value == "volume") {
				APF.scaling = NULL;
			}
			else {
				for (i = 0; i < dimension; i++) {
					APF.scaling[i] = APF.Bound_N[i] - APF.Bound_I[i];
					if (APF.scaling[i] == 0) {
						printf("WARNING! COORDINATE %d HAS THE SAME VALUE FOR IDEAL AND BOUND FOR NADIR POINT", i);
						APF.scaling[i] = DELTA;
						warning = true;
					}
				}
			}
			printf("DONE.");
			return true;
		}
};







