vector<point> tangentsPointPolygon(const vector<point> & P, const vector<vector<point>> & Ps, const point & p){
	int n = P.size(), m = Ps[0].size(), k = Ps[1].size();
	
	int lk = m; if(Ps[0][m - 1] == Ps[1][0]) lk--; 

	auto tang = [&](int l, int r, ld w, int kl) -> int {
		int res = min(l, r);
	        while(l <= r){
			int m = (l + r) / 2;
			ld a = (P[(m + kl) % n] - p).cross(P[(m + 1 + kl) % n] - p) * w, b = (P[(m + kl) % n] - p).cross(P[(m - 1 + n + kl) % n] - p) * w;
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
			else res = m, l = m + 1;
		}
	        return res;
	    };
    
	point left = p, rigth = p;

	int t1 = bs(0, m - 1, Ps[0], 1), t2 = bs(0, k - 1, Ps[1], -1);
	
	auto u1 = tang(0, t1, -1, 0), u2 = tang(0, t2, -1, lk);
	auto v1 = tang(t1, m - 1, 1, 0), v2 = tang(t2, k - 1, 1, lk);
	
	if(leq((P[u1] - p).cross(P[(u1 - 1 + n) % n] - p), 0) && leq((P[u1] - p).cross(P[(u1 + 1) % n] - p), 0)) left = P[u1];
	if(leq((P[(lk + u2) % n] - p).cross(P[(lk + u2 - 1 + n) % n] - p), 0) && leq((P[(lk + u2) % n] - p).cross(P[(lk + u2 + 1) % n] - p), 0)) left = P[(lk + u2) % n];
	
	if(geq((P[v1] - p).cross(P[(v1 - 1 + n) % n] - p), 0) && geq((P[v1] - p).cross(P[(v1 + 1) % n] - p), 0)) rigth = P[v1];
	if(geq((P[(lk + v2) % n] - p).cross(P[(lk - 1 + n + v2) % n] - p), 0) && geq((P[(lk + v2) % n] - p).cross(P[(lk + 1 + v2) % n] - p), 0)) rigth = P[(lk + v2) % n];
    
	return {left, rigth};
}
