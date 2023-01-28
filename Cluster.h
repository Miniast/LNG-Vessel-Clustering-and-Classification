#pragma GCC optimize(2)
#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#include<cstdio>
#include<vector>
#include<ctime>
#include<unordered_map>
#include<cmath>
#include "csv.h"
#define R register
using namespace std;
typedef long long ll;
const int MAXN = 5e6+5;
const int MAXM = 500005;
const int MAXK = 1005;
const int end_of_file = 998244353;
const double inf = 1000.0;

struct node{
	int mmsi,unix_time,draft;
}vessel[MAXN];

struct nodes{
	int number;
    double longtitude,latitude;
	vector <int> node_list;
}batch[MAXM];

struct cluster{
	double avg_longti,avg_lati,sum_longti,sum_lati;
	int number,type;
	vector <int> batch_list;
}clst[MAXK];

unordered_map <ll,int> mp;
int n,m,times;
int K;
int mylist[MAXM],ord[MAXM];
