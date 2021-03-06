#include <iostream> 
#include<string>
#include<vector>
#include <algorithm>
using namespace std;


int rank_[1000001];
int tmp[1000001];


//文字列Sの接尾辞配列を構築
void construct_sa(string S, vector<int>&sa) {
	int n = S.length();
	
	auto compare = [n](int k) {
		return [n,k](int i, int j) {
			if (rank_[i] != rank_[j]) {
				return rank_[i] < rank_[j];
			}
			else {
				int ri = i + k <= n ? rank_[i + k] : -1;
				int rj = j + k <= n ? rank_[j + k] : -1;

				return ri < rj;
			}
		};
	};
	for (int i = 0; i <= n; i++) {//0~nまで
		sa[i] = i;
		rank_[i] = i < n ? S[i] : -1;//全てのk文字の部分文字列をソートしたときにS[i,k]が何番目に小さいか表す
		//iがnのときrankを-1で埋める
	}
	//k文字についてソートされているところから,2k文字でソートする
	for (int k = 1; k <= n; k *= 2) {
		sort(sa.begin(), sa.end(), compare(k));

		//いったんtmpに次にランクを計算し,それからrankに移す
		tmp[sa[0]] = 0;//空文字列のrankは常に0
		for (int i = 1; i <= n; i++) {

			tmp[sa[i]] = tmp[sa[i - 1]] + (compare(k)(sa[i - 1], sa[i]) ? 1 : 0);
		}
		for (int i = 0; i <= n; i++) {
			rank_[i] = tmp[i];
		}
	}
	

}
vector<string>rebuild(vector<int>& sa, string T) {
	vector<string>res;
	for (int i = 0; i < sa.size(); i++) {
		string T_sub = T.substr(sa[i], T.length());
		res.push_back(T_sub);
	}
	return res;

}

bool contain(string S, vector<int>& sa, string T) {
	int a = 0;
	int b = S.length();
	while (b - a>1) {
		int c = (a + b) / 2;
		string S_sub = S.substr(sa[c], T.length());
		if (S_sub.compare(T) < 0) {
			a = c;
		}
		else {
			b = c;
		}
		
	}
	string S_sub2 = S.substr(sa[b], T.length());
	return S_sub2.compare(T) == 0;
}

int main() {
	//abracadabra
	string T;
	cin >> T;
	int n = T.length();
	vector<int>sa(n+1,0);
	construct_sa(T, sa);
	/*
	auto res = rebuild(sa, T);
	for (int i = 0; i < res.size(); i++) {
		cout << res[i] << endl;
	}
	*/
	
	
	
	int Q;
	cin >> Q;
	for (int i = 0; i < Q; i++) {
		string P;
		cin >> P;
		int ans = contain(T, sa, P);
		cout << ans << endl;

	}
	
	
	
	return 0;
}