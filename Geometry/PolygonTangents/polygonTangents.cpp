vector<point> tangentsPointPolygon(const vector<point> & P, const vector<vector<point>> & Ps, const point & p, const vector<vector<point>> & T, const vector<vector<int>> & Ls, const vector<point> & B){
	for(int i = 0; i < T.size(); i++){
		if(pointInPolygon(T[i], p)){
			return {T[i][1], T[i][0]}; 
		}
	}

	int n = P.size(), m = Ps[0].size(), k = Ps[1].size();

	auto tang = [&](int l, int r, const vector<point> & A, ld w) -> int {
		int res = l;
		int y = A.size();
	        while(l <= r){
	        	int m = (l + r) / 2;
				ld a = (A[m] - p).cross(A[(m + 1) % y] - p) * w, b = (A[m] - p).cross(A[(m - 1 + y) % y] - p) * w;
				if(geq(a, 0) && geq(b, 0)) return m;
	            	if(geq(a, 0)) r = m - 1, res = m;
	            	else l = m + 1;
		}
	        return res;
	};

	auto bs = [&](int l, int r, const vector<point> & A, ld w) -> int {
	        int res = l;
	        ld w1 = p.x * w;
	        while(l <= r){
	            	int m = (l + r) / 2;
				if(ge(A[m].x * w, w1)) r = m - 1;
				else res = max(res, m), l = m + 1;
	        }
        return res;
    	};

	point left = p, rigth = p;
	point f1 = Ps[0][0], f2 = Ps[0][m - 1];
	if(neq((P[Ls[0][1]] - P[Ls[0][0]]).cross(P[Ls[1][1]] - P[Ls[1][0]]), 0)) f1 = intersectLines(P[Ls[0][0]], P[Ls[0][1]] - P[Ls[0][0]], P[Ls[1][0]], P[Ls[1][1]] - P[Ls[1][0]]);
	if(neq((P[Ls[2][1]] - P[Ls[2][0]]).cross(P[Ls[3][1]] - P[Ls[3][0]]), 0)) f2 = intersectLines(P[Ls[2][0]], P[Ls[2][1]] - P[Ls[2][0]], P[Ls[3][0]], P[Ls[3][1]] - P[Ls[3][0]]);
	
	if(geq((P[Ls[0][1]] - P[Ls[0][0]]).cross(p - P[Ls[0][0]]), 0) && geq((P[Ls[1][1]]- P[Ls[1][0]]).cross(p - P[Ls[1][0]]), 0) && Ls[0][0] != Ls[0][1]){
		left = Ps[1][tang(0, k - 1, Ps[1], -1)];
		rigth = Ps[0][tang(0, m - 1, Ps[0], 1)];
	}else if(geq((P[Ls[2][1]] - P[Ls[2][0]]).cross(p - P[Ls[2][0]]), 0) && geq((P[Ls[3][1]]- P[Ls[3][0]]).cross(p - P[Ls[3][0]]), 0) && Ls[2][0] != Ls[2][1]){
		left = Ps[0][tang(0, m - 1, Ps[0], -1)];
		rigth = Ps[1][tang(0, k - 1, Ps[1], 1)];
	}else if(le((f2 - f1).cross(p - f1), 0)){
		int t = bs(0, m - 1, Ps[0], 1);
		if(leq((Ps[1][k - 1] - p).cross(Ps[0][0] - p), 0) && leq((Ps[1][k - 1] - p).cross(Ps[1][k - 2] - p), 0) && leq((Ps[1][k - 1] - p).cross(Ps[0][1] - p), 0)) left = Ps[1][k - 1];
		else left = Ps[0][tang(0, min(t, m - 2), Ps[0], -1)];
		if(geq((Ps[1][0] - p).cross(Ps[0][m - 1] - p), 0) && geq((Ps[1][0] - p).cross(Ps[1][1] - p), 0) && geq((Ps[1][0] - p).cross(Ps[0][m - 2] - p), 0)) rigth = Ps[1][0];
		else rigth = Ps[0][tang(min(t + 1, m - 1), m - 1, Ps[0], 1)];
	}else if(ge((f2 - f1).cross(p - f1), 0)){
		int t = bs(0, k - 1, Ps[1], -1);
		if(leq((Ps[0][m - 1] - p).cross(Ps[1][0] - p), 0) && leq((Ps[0][m - 1] - p).cross(Ps[0][m - 2] - p), 0) && leq((Ps[0][m - 1] - p).cross(Ps[1][1] - p), 0)) left = Ps[0][m - 1];
		else left = Ps[1][tang(0, min(t, k - 2), Ps[1], -1)];
		if(geq((Ps[0][0] - p).cross(Ps[1][k - 1] - p), 0) && geq((Ps[0][0] - p).cross(Ps[1][k - 2] - p), 0) && geq((Ps[0][0] - p).cross(Ps[0][1] - p), 0)) rigth = Ps[0][0];
		else rigth = Ps[1][tang(min(t + 1, k - 1), k - 1, Ps[1], 1)];
	}else{
		int i = lower_bound(all(Ps[0]), p) - Ps[0].begin();
		int j = lower_bound(all(B), p) - B.begin();
		if(pointInSegment(Ps[0][i], Ps[0][i - 1], p)) left = Ps[0][i - 1], rigth = Ps[0][i];
		if(pointInSegment(B[j], B[j - 1], p)) left = B[j], rigth = B[j - 1];
	}
	return {left, rigth};
}

pair<vector<vector<point>>, vector<vector<int>>> precTangents(vector<point> P){
	int n = P.size();
	auto R = twoSidesCH(P);
	int m = R[0].size(), k = R[1].size();
	bool f1 = false, f2 = false;
	int il1 = 0, il2 = 0, sl1 = 0, sl2 = 0, ir1 = 0, ir2 = 0, sr1 = 0, sr2 = 0;
	il1 = 1; il2 = 0; sl1 = 0; sl2 = n - 1;
	ir1 = m - 1; ir2 = m - 2; sr1 = m % n; sr2 = m - 1;
	if(R[0][m - 1] != R[1][0]) sr1 = (sr1 + 1) % n, sr2 = (sr2 + 1) % n;
	if(R[0][0] != R[1][k - 1]) sl1 = n - 1, sl2 = n - 2;
	point a(-INF2, -1000), b(-INF2, 1000);
	vector<vector<point>> Ts;
	if(il2 != sl1){
		if(eq((P[il2] - P[il1]).cross(P[sl2] - P[sl1]), 0)){
			auto it = intersectLines(P[il1], P[il2] - P[il1], a, b - a), at = intersectLines(P[sl1], P[sl2] - P[sl1], a, b - a);
			Ts.pb({P[il2], P[sl1], at, it});
			il1 = il2 = sl1 = sl2 = -1;
		}else{
			auto it = intersectLines(P[il1], P[il2] - P[il1], P[sl1], P[sl2] - P[sl1]);
			if(ge(it.x, P[il1].x)){
				auto it = intersectLines(P[il1], P[il2] - P[il1], a, b - a), at = intersectLines(P[sl1], P[sl2] - P[sl1], a, b - a);
				Ts.pb({P[il2], P[sl1], at, it});
				il1 = il2 = sl1 = sl2 = -1;
			}else{
				Ts.pb({P[il2], P[sl1], it});
			}
		}
	}
	a = point(INF2, -1000); b = point(INF2, 1000);
	if(sr2 != ir1){
		if(eq((P[ir2] - P[ir1]).cross(P[sr2] - P[sr1]), 0)){
			auto it = intersectLines(P[ir1], P[ir2] - P[ir1], a, b - a), at = intersectLines(P[sr1], P[sr2] - P[sr1], a, b - a);
			Ts.pb({P[sr2], P[ir1], it, at});
			ir1 = ir2 = sr1 = sr2 = -1;
		}else{
			auto it = intersectLines(P[ir1], P[ir2] - P[ir1], P[sr1], P[sr2] - P[sr1]);
			if(le(it.x, P[ir1].x)){
				auto it = intersectLines(P[ir1], P[ir2] - P[ir1], a, b - a), at = intersectLines(P[sr1], P[sr2] - P[sr1], a, b - a);
				Ts.pb({P[sr2], P[ir1], it, at});
				ir1 = ir2 = sr1 = sr2 = -1;
			}else{
				Ts.pb({P[sr2], P[ir1], it});
			}
		}
	}
	return {Ts, {{il1, il2}, {sl1, sl2}, {ir1, ir2}, {sr1, sr2}}};
}
