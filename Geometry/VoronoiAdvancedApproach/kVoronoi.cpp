vector<vector<vector<point>>> kVoronoi(vector<point> P, int k){

	int n = P.size();

	auto subVor = [&](vector<point> A, int act, vector<int> neigh) -> vector<point> {
		int m = neigh.size();
		for(int i = 0; i < m; i++){
			if(act == neigh[i]) continue;
			A = cutPolygon(A, (P[act] + P[neigh[i]]) / 2, (P[act] - P[neigh[i]]).perp());
		}
		return A;
	};

	auto identify = [&](vector<point> A, int act) -> vector<int> {
		vector<int> neigh;
		int m = A.size();
		for(int i = 0; i < n; i++){
			if(i == act) continue;
			for(int j = 0; j < m; j++){
				if(eq((P[act] - (A[j] - A[(j + 1) % m])).length(), (P[i] - (A[j] - A[(j + 1) % m])).length())) neigh.pb(i);
			}
		}
		return neigh;
	};

	auto inicialice = [&]() -> vector<point> {
		return {point(-Inf, -Inf), point(Inf, -Inf), point(Inf, Inf), point(-Inf, Inf)};
	};

	vector<vector<vector<point>>> ans(n), res(n);
	vector<vector<vector<int>>> neighs(n);

	for(int i = 0; i < n; i++){
		vector<int> u;
		for(int j = 0; j < n; j++) if(i != j) u.pb(j);
		ans[i].pb(subVor(inicialice(), i, u));
	}

	for(int i = 1; i < k; i++){

		for(int j = 0; j < n; j++){
			if(i == 0) continue;
			for(int l = 0; l < ans[j].size(); l++){
				neighs[j].pb(identify(ans[j][l], j));
			}
		}

		for(int j = 0; j < n; j++){
			for(int l = 0; l < ans[j].size(); l++){
				for(int t = 0; t < neighs[j][l].size(); t++){
					res[t].pb(subVor(ans[j][l], neighs[j][l][t], neighs[j][l]));
				}
			}
		}

		ans = res;

		res.clear();
		neighs.clear();

	}

	return ans;

}
