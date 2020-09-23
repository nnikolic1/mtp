#include <vector>
#include <iostream>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <chrono>
#include <ctime>
#include <cmath>

unsigned timeNow() {
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

using namespace std; 

const unsigned n = 2000;
const long double t = 4000;
unsigned m;
// Ako se svaka moguća grana dodaje u graf sa verovatnoćom p,
// svaka neuređena trojka čvorova će činiti trougao sa verovatnoćom p^3,
// pa je očekivani broj trouglova (n*(n-1)*(n-2)/6)*p^3
const long double p = cbrt(6*t/n/(n-1)/(n-2));
bool adj[n][n];
vector< tuple<int, int, int> > triangles;
unsigned lb = 0;
vector<bool> BnBopt;

bool triangles_intersect(const tuple<int, int, int>& i, const tuple<int, int, int>& j) {
	return (get<0>(i)==get<0>(j)) || (get<0>(i)==get<1>(j)) || (get<0>(i)==get<2>(j))
		|| (get<1>(i)==get<0>(j)) || (get<1>(i)==get<1>(j)) || (get<1>(i)==get<2>(j))
		|| (get<2>(i)==get<0>(j)) || (get<2>(i)==get<1>(j)) || (get<2>(i)==get<2>(j));
}

vector< tuple<int, int, int> > LBbeta(vector< tuple<int, int, int> > x) {
	vector< tuple<int, int, int> > y;
	unsigned i,m1,minI,minB;
	tuple<int, int, int> t;
	while (!x.empty()) {
		m1 = x.size();
		// alpha[i] je broj trouglova iz liste kojima pripada čvor sa indeksom i
		// beta[i] je težina x[i] tj. zbir alpha[i] njegova tri čvora
		vector<unsigned> alpha(n,0);
		vector<unsigned> beta(m1);
		for (i=0;i<m1;i++) {
			alpha[get<0>(x[i])]++;
			alpha[get<1>(x[i])]++;
			alpha[get<2>(x[i])]++;
		}
		// Određujemo trougao u listi sa minimalnom težinom, njegov indeks je minI
		beta[0] = alpha[get<0>(x[0])]+alpha[get<1>(x[0])]+alpha[get<2>(x[0])];
		minB = beta[0];
		minI = 0;
		for (i=1;i<m1;i++) {
			beta[i] = alpha[get<0>(x[i])]+alpha[get<1>(x[i])]+alpha[get<2>(x[i])];
			if (beta[i] < minB) {
				minB = beta[i];
				minI = i;
			}
		}
		// Uzimamo trougao u rešenje i iz liste brišemo njega i njegove susede
		t = x[minI];
		y.push_back(t);
		x.erase(remove_if(x.begin(),x.end(),[t](tuple<int, int, int> a){return triangles_intersect(t,a);}),x.end());
	}
	return y;
}

// Implementacije algoritama iz literature

vector< tuple<int, int, int> > LB1(vector< tuple<int, int, int> > x) {
	vector< tuple<int, int, int> > y;
	unsigned i,j,m1,minI,minL;
	tuple<int, int, int> t;
	while (!x.empty()) {
		m1 = x.size();
		vector<unsigned> lambda(m1,0);
		for (i=0;i<m1;i++)
			for (j=i+1;j<m1;j++)
				if (triangles_intersect(x[i],x[j])) {
					lambda[i]++;
					lambda[j]++;
				}
		minI = 0;
		minL = lambda[0];
		for (i=1;i<m1;i++)
			if (lambda[i] < minL) {
				minI = i;
				minL = lambda[i];
			}
		t = x[minI];
		y.push_back(t);
		x.erase(remove_if(x.begin(),x.end(),[t](tuple<int, int, int> a){return triangles_intersect(t,a);}),x.end());
	}
	return y;
}

vector<bool> knapsackEqualValues(vector< pair<int,int> > x, unsigned max) {
	vector<bool> r(m,0);
	while (1) {
		if (x.size() == 0)
			return r;
		unsigned min = (*min_element(x.begin(),x.end(),[](pair<int,int> i,pair<int,int> j){return i.second<j.second;})).second;
		if (min > max)
			return r;
		auto pivot = x[rand()%x.size()];
		vector < pair<int,int> > xl,xd;
		vector<int> li;
		unsigned lsum = 0;
		for (auto s: x)
			if (s.second < pivot.second || (s.second == pivot.second && s.first <= pivot.first)) {
				lsum += s.second;
				xl.push_back(s);
				li.push_back(s.first);
			}
			else
				xd.push_back(s);
		if (lsum <= max) {
			x = xd;
			max -= lsum;
			for (auto s: li)
				r[s] = 1;
		}
		else
			x = xl;
	}
}

/*vector<bool> LB2() {
	vector<bool> x(m);
	bool solved;
	unsigned i,j,k;
	vector<bool> removed(m,0);
	do {
		x = vector<bool>(m,0);
		auto max = 3*count_if(removed.begin(),removed.end(),[](bool b){return !b;});
		vector<unsigned> alpha(n,0);
		for (i=0;i<m;i++)
			if (!removed[i]) {
				alpha[get<0>(triangles[i])]++;
				alpha[get<1>(triangles[i])]++;
				alpha[get<2>(triangles[i])]++;
			}
		vector< pair<int,int> > zip(m);
		for (i=0;i<m;i++)
			if (!removed[i])
				zip[i] = make_pair(i,alpha[get<0>(triangles[i])]+alpha[get<1>(triangles[i])]+alpha[get<2>(triangles[i])]);
			else
				zip[i] = make_pair(i,max+1);
		x = knapsackEqualValues(zip,max);
		solved = 1;
		bool br = 0;
		unsigned c;
		for (i=0;i<m && !br;i++)
			for (j=i+1;j<m && !br;j++)
				if (x[i] && x[j] && triangles_intersect(i,j)) {
					auto betaI = alpha[get<0>(triangles[i])]+alpha[get<1>(triangles[i])]+alpha[get<2>(triangles[i])];
					auto betaJ = alpha[get<0>(triangles[j])]+alpha[get<1>(triangles[j])]+alpha[get<2>(triangles[j])];
					c = (betaJ < betaI) ? j : i;
					for (k=0;k<m;k++)
						if (k != c && triangles_intersect(c,k))
							removed[k] = 1;
					solved = 0;
					br = 1;
				}
	} while (!solved);
	return x;
}

unsigned UB2(const vector<bool>& removed) {
	unsigned i;
	auto max = 3*count_if(removed.begin(),removed.end(),[](bool b){return !b;});
	vector<unsigned> alpha(n,0);
	for (i=0;i<m;i++)
		if (!removed[i]) {
			alpha[get<0>(triangles[i])]++;
			alpha[get<1>(triangles[i])]++;
			alpha[get<2>(triangles[i])]++;
		}
	vector< pair<int,int> > zip(m);
	for (i=0;i<m;i++)
		if (!removed[i])
			zip[i] = make_pair(i,alpha[get<0>(triangles[i])]+alpha[get<1>(triangles[i])]+alpha[get<2>(triangles[i])]);
		else
			zip[i] = make_pair(i,max+1);
	auto kn = knapsackEqualValues(zip,max);
	return accumulate(kn.begin(),kn.end(),0);
}

unsigned UB1(const vector<bool>& removed) {
	unsigned i,incl = 0;
	vector<bool> included(n,0);
	for (i=0;i<m;i++) {
		if (!removed[i]) {
			if (!included[get<0>(triangles[i])]) {
				included[get<0>(triangles[i])] = 1;
				incl++;
			}
			if (!included[get<1>(triangles[i])]) {
				included[get<1>(triangles[i])] = 1;
				incl++;
			}
			if (!included[get<2>(triangles[i])]) {
				included[get<2>(triangles[i])] = 1;
				incl++;
			}
		}
	}
	return incl/3;
}

unsigned UB(const vector<bool>& removed) {
	auto ub1 = UB1(removed);
	auto ub2 = UB2(removed);
	return (ub1 < ub2) ? ub1 : ub2;
}

vector<bool> GS() {
	vector<bool> x(m,0);
	unsigned i,j,ub=0;
	vector<bool> restricted(m,0);
	for (i=0;i<m;i++) {
		auto restrictedL = restricted;
		restrictedL[i] = 1;
		auto ubL = ub + UB(restrictedL);
		bool solvedR = 1;
		auto restrictedR = restricted;
		for (j=0;j<m;j++)
			if (triangles_intersect(i,j)) {
				if (j < i && x[j]) {
					solvedR = 0;
					break;
				}
				else
					restrictedR[j] = 1;
			}
		unsigned ubR = 0;
		if (solvedR)
			ubR = ub + 1 + UB(restrictedR);
		if (ubL < ubR) {
			ub++;
			restricted = restrictedR;
			x[i] = 1;
		}
		else
			restricted = restrictedL;
	}
	return x;
}

vector<bool> LB() {
	auto lb1 = LB1();
	auto lb1T = accumulate(lb1.begin(),lb1.end(),0);
	auto lb2 = LB2();
	auto lb2T = accumulate(lb2.begin(),lb2.end(),0);
	auto lb3 = GS();
	auto lb3T = accumulate(lb3.begin(),lb3.end(),0);
	return (lb1T > lb2T) ? ((lb1T > lb3T) ? lb1 : lb3) : ((lb2T > lb3T) ? lb2 : lb3);
}

void BnB(const vector<bool>& restricted, vector<bool> x, unsigned i = 0, unsigned ub = 0) {
	if (ub > lb) {
		cout << "lb = " << ub << endl;
		lb = ub;
		BnBopt = x;
	}
	if (i >= m)
		return;
	unsigned j;
	auto restrictedL = restricted;
	restrictedL[i] = 1;
	auto ubL = ub + UB(restrictedL);
	bool solvedR = 1;
	auto restrictedR = restricted;
	for (j=0;j<m;j++)
		if (triangles_intersect(i,j)) {
			if (j < i && x[j]) {
				solvedR = 0;
				break;
			}
			else
				restrictedR[j] = 1;
		}
	unsigned ubR = 0;
	if (solvedR)
		ubR = ub + 1 + UB(restrictedR);
	if (ubL > lb && ubR > lb) {
		if (ubL > ubR) {
			BnB(restrictedL,x,i+1,ub);
			x[i] = 1;
			BnB(restrictedR,x,i+1,ub+1);
		}
		else {
			x[i] = 1;
			BnB(restrictedR,x,i+1,ub+1);
			x[i] = 0;
			BnB(restrictedL,x,i+1,ub);
		}
	}
	else if (ubL > lb) {
		BnB(restrictedL,x,i+1,ub);
	}
	else if (ubR > lb) {
		x[i] = 1;
		BnB(restrictedR,x,i+1,ub+1);
	}
}

vector<bool> randomVector() {
	vector<bool> s(m,0);
	unsigned i;
	for (i=0;i<m;i++)
		if ((long double)rand()/RAND_MAX < 0.03)
			s[i] = 1;
	return s;
}

vector<bool> randomSolution() {
	unsigned i,j;
	bool intersect;
	vector<bool> s(m,0);
	vector<unsigned> ind(m);
	for (i=0;i<m;i++)
		ind[i] = i;
	random_shuffle(ind.begin(),ind.end());
	for (i=0;i<m;i++) {
		intersect = 0;
		for (j=0;j<i;j++)
			if (s[ind[j]] && triangles_intersect(ind[i],ind[j])) {
				intersect = 1;
				break;
			}
		if (!intersect)
			s[ind[i]] = 1;
	}
	return s;
}

unsigned fitness(vector<bool> x) {
	unsigned i,j;
	for (i=0;i<m;i++)
		if (x[i])
			for (j=i+1;j<m;j++)
				if (x[j] && triangles_intersect(i,j))
					return 0;
	return accumulate(x.begin(),x.end(),0);
}

vector< pair< vector<bool>,unsigned > > generatePopulation() {
	vector< pair< vector<bool>,unsigned > > population(100);
	BnBopt = vector<bool>(m,0);
	unsigned i;
	for (i=0;i<population.size();i++) {
		pair< vector<bool>,unsigned > x;
		auto temp = randomVector();
		x.first = temp;
		x.second = fitness(temp);
		population[i] = x;
	}
	sort(population.begin(),population.end(),[](pair< vector<bool>,unsigned > x,pair< vector<bool>,unsigned > y){return x.second>y.second;});
	if (population[0].second > lb) {
		lb = population[0].second;
		BnBopt = population[0].first;
	}
	return population;
 }

void GA(vector< pair< vector<bool>,unsigned > >& population) {
	unsigned popSize = population.size();
	unsigned withoutImprovement = 0;
	while (withoutImprovement < 100) {
		unsigned s,r,j;
		auto p1 = population[rand()%(popSize/2)].first;
		auto p2 = population[(popSize/2)+rand()%(popSize/2)].first;
		vector<bool> c1 = vector<bool>(m);
		vector<bool> c2 = vector<bool>(m);
		r = rand()%m;
		s = rand()%m;
		if (s < r) {
			auto temp = s;
			s = r;
			r = temp;
		}
		for (j=0;j<r;j++) {
			c1[j] = p1[j];
			c2[j] = p2[j];
		}
		for (;j<s;j++) {
			c1[j] = p2[j];
			c2[j] = p1[j];
		}
		for (;j<m;j++) {
			c1[j] = p1[j];
			c2[j] = p2[j];
		}
		auto rm = ((long double)rand()/RAND_MAX);
		if (rm < 0.05) {
			r = rand()%m;
			c1[r] = !c1[r];
		}
		rm = ((long double)rand()/RAND_MAX);
		if (rm < 0.05) {
			r = rand()%m;
			c2[r] = !c2[r];
		}
		auto c1p = make_pair(c1,fitness(c1));
		auto c2p = make_pair(c2,fitness(c2));
		population[popSize-2] = c1p;
		population[popSize-1] = c2p;
		sort(population.begin(),population.end(),[](pair< vector<bool>,unsigned > x,pair< vector<bool>,unsigned > y){return x.second>y.second;});
		withoutImprovement++;
		if (population[0].second > lb) {
			lb = population[0].second;
			BnBopt = population[0].first;
			withoutImprovement = 0;
		}
	}
}*/

int main() {
	unsigned i,j,k,t1,t2;
	srand(timeNow());
	// Slučajno generišemo graf
	// Svaka moguća grana se dodaje sa verovatnoćom p
	for (i=0;i<n;i++)
		for (j=i+1;j<n;j++)
			if (((long double)rand()/RAND_MAX) < p) {
				adj[i][j] = 1;
				adj[j][i] = 1;
			}
			else {
				adj[i][j] = 0;
				adj[j][i] = 0;
			}
	// Nalazimo sve trouglove u grafu
	t1 = timeNow();
	for (i=0;i<n;i++)
		for (j=i+1;j<n;j++)
			if (adj[i][j])
				for (k=j+1;k<n;k++)
					if (adj[j][k] && adj[k][i])
						triangles.push_back(make_tuple(i,j,k));
	t2 = timeNow();
	m = triangles.size();
	vector< tuple<int, int, int> > temp;
	cout << m << " trouglova za " << t2-t1 << "ms\n";
	t1 = timeNow();
	temp = LB1(triangles);
	t2 = timeNow();
	cout << "LB1:    " << temp.size() << " za " << t2-t1 << "ms\n";
	t1 = timeNow();
	temp = LBbeta(triangles);
	t2 = timeNow();
	cout << "LBbeta: " << temp.size() << " za " << t2-t1 << "ms\n";
	/*t1 = timeNow();
	LB4();
	t2 = timeNow();
	cout << "LB4: " << lb << " za " << t2-t1 << "ms\n";*/
	/*t1 = timeNow();
	temp = GS();
	t2 = timeNow();
	cout << "LB3: " << accumulate(temp.begin(),temp.end(),0) << " za " << t2-t1 << "ms\n";
	t1 = timeNow();
	auto pop = generatePopulation();
	t2 = timeNow();
	cout << "pop: " << lb << " za " << t2-t1 << "ms\n";
	t1 = timeNow();
	GA(pop);
	t2 = timeNow();
	cout << "GA: " << lb << " za " << t2-t1 << "ms\n";*/
	return 0;
}