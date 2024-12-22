#include <bits/stdc++.h>
using namespace std;

// Holi c:

#define ll long long int
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()

const int Inf = 1e3;
const ll mod = 1e9+7;
const ll INF = 1e18;

using ld = long double;
const ld eps = 1e-9, inf = numeric_limits<ld>::max(), pi = acos(-1);
bool geq(ld a, ld b){return a-b >= -eps;}     //a >= b
bool leq(ld a, ld b){return b-a >= -eps;}     //a <= b
bool ge(ld a, ld b){return a-b > eps;}        //a > b
bool le(ld a, ld b){return b-a > eps;}        //a < b
bool eq(ld a, ld b){return abs(a-b) <= eps;}  //a == b
bool neq(ld a, ld b){return abs(a-b) > eps;}  //a != b

struct point3{
    ld x, y, z;
    point3(): x(0), y(0), z(0){}
    point3(ld x, ld y, ld z): x(x), y(y), z(z){}

    point3 operator+(const point3 & p) const{return point3(x + p.x, y + p.y, z + p.z);}
    point3 operator-(const point3 & p) const{return point3(x - p.x, y - p.y, z - p.z);}
    point3 operator*(const ld & k) const{return point3(x * k, y * k, z * k);}
    point3 operator/(const ld & k) const{return point3(x / k, y / k, z / k);}

    point3 operator+=(const point3 & p){*this = *this + p; return *this;}
    point3 operator-=(const point3 & p){*this = *this - p; return *this;}
    point3 operator*=(const ld & p){*this = *this * p; return *this;}
    point3 operator/=(const ld & p){*this = *this / p; return *this;}

    ld dot(const point3 & p) const{return x * p.x + y * p.y + z * p.z;}
    point3 cross(const point3 & p) const{return {y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x};}
    ld norm() const{return x * x + y * y + z * z;}
    ld length() const{return sqrt(x * x + y * y + z * z);}
    point3 unit() const{return (*this) / length();}
    ld angle(const point3 & p) const{ld a = atan2l(((*this).cross(p)).norm(), (*this).dot(p)); if(le(a, 0)) a += pi * 2; return a;}

    bool operator==(const point3 & p) const{return eq(x, p.x) && eq(y, p.y) && eq(z, p.z);}
    bool operator!=(const point3 & p) const{return !(*this == p);}
    bool operator<(const point3 & p) const{return le(x, p.x) || (eq(x, p.x) && le(y, p.y)) || (eq(x, p.x) && eq(y, p.y) && le(z, p.z));}
    bool operator>(const point3 & p) const{return ge(x, p.x) || (eq(x, p.x) && ge(y, p.y)) || (eq(x, p.x) && eq(y, p.y) && ge(z, p.z));}
    bool operator<=(const point3 & p) const{return !(*this > p);}
    bool operator>=(const point3 & p) const{return !(*this < p);}
};

istream &operator>>(istream &is, point3 & p){return is >> p.x >> p.y >> p.z;}
ostream &operator<<(ostream &os, const point3 & p){return os << "(" << p.x << ", " << p.y << " , " << p.z << ")";}

struct plane{
    point3 n; ld d;
    plane(): n(0, 0, 0), d(0){}
    plane(point3 n, ld d): n(n), d(d){}
    plane(point3 n, point3 p): n(n), d(n.dot(p)){}
    plane(point3 p1, point3 p2, point3 p3): plane((p2 - p1).cross(p3 - p1), p1.dot((p2 - p1).cross(p3 - p1))){}

    ld side(const point3 & p) const{return ((*this).n).dot(p) - (*this).d;}
    ld dist(const point3 & p) const{return abs((*this).side(p)) / ((*this).n).length();}
    plane translate(const point3 & p) const{return {(*this).n, (*this).d + ((*this).n).dot(p)};}
    plane shift(const ld & t) const{return {(*this).n, (*this).d + t * ((*this).n).length()};}
    point3 proj(const point3 & p) const{return p - (*this).n * (*this).side(p) / ((*this).n).length();}
    point3 refl(const point3 & p) const{return p - (*this).n * (*this).side(p) / ((*this).n).length() * 2;}
    //Minium angle between two planes
    ld anglep(const plane & p) const{return ((*this).n).angle(p.n);}
    //Angle between a plane and line is pi / 2 - angle(between two vector directionals) ans parallel and perpendicular is exactly of plane to plane

    bool operator|(const plane & p) const{const point3 r = ((*this).n).cross(p.n); return (r.x == 0 && r.y == 0 && r.z == 0);} //Parallel
    bool operator/(const plane & p) const{return ((*this).n).dot(p.n) == 0;} //Perpendicular
};

struct line{
    point3 p, q;
    line(): p(0, 0, 0), q(0, 0, 0){}
    line(point3 p, point3 q): p(p), q(q){}

    ld dist(const point3 & p) const{return sqrt((((*this).q).cross(((*this).p) - p)).norm() / ((*this).q).norm());}
    point3 proj(const point3 & p) const{return (*this).p + ((*this).q) * (((*this).q).dot(p - (*this).p)) / ((*this).q).norm();}
    point3 refl(const point3 & p) const{return (*this).proj(p) * 2 - p;}
    //assuming intersection exist and p.n and d not collinear (p.n).dot(d) != 0
    point3 inter(const plane & p) const{return (*this).p - ((*this).q) * (p.side((*this).p)) / (p.n).dot((*this).q);}
    //Minium angle between two lines
    ld angle(const line & l) const{return ((*this).q).angle(l.q);}
    //Angle between a line and plane is pi / 2 - angle(between two vector directionals) ans parallel and perpendicular is exactly of line to line

    bool operator|(const line & l) const{const point3 r = ((*this).q).cross(l.q); return (r.x == 0 && r.y == 0 && r.z == 0);} //Parallel
    bool operator/(const line & l) const{return ((*this).q).dot(l.q) == 0; } //Parallel
};

struct point{
	ld x, y;
	point(): x(0), y(0){}
	point(ld x, ld y): x(x), y(y){}

	point operator+(const point & p) const{return point(x + p.x, y + p.y);}
	point operator-(const point & p) const{return point(x - p.x, y - p.y);}
	point operator*(const ld & k) const{return point(x * k, y * k);}
	point operator/(const ld & k) const{return point(x / k, y / k);}

	point operator+=(const point & p){*this = *this + p; return *this;}
	point operator-=(const point & p){*this = *this - p; return *this;}
	point operator*=(const ld & p){*this = *this * p; return *this;}
	point operator/=(const ld & p){*this = *this / p; return *this;}
	point perp() const{return point(-y, x);}
	ld ang() const{
		ld a = atan2l(y, x); a += le(a, 0) ? 2*pi : 0; return a;
	}
	ld dot(const point & p) const{return x * p.x + y * p.y;}
	ld cross(const point & p) const{return x * p.y - y * p.x;}
	ld norm() const{return x * x + y * y;}
	ld length() const{return sqrtl(x * x + y * y);}
	point unit() const{return (*this) / length();}

	bool operator==(const point & p) const{return eq(x, p.x) && eq(y, p.y);}
	bool operator!=(const point & p) const{return !(*this == p);}
	bool operator<(const point & p) const{return le(x, p.x) || (eq(x, p.x) && le(y, p.y));}
	bool operator>(const point & p) const{return ge(x, p.x) || (eq(x, p.x) && ge(y, p.y));}
	bool half(const point & p) const{return le(p.cross(*this), 0) || (eq(p.cross(*this), 0) && le(p.dot(*this), 0));}
};

void polarSort(vector<point> & P, const point & o, const point & v){
	//sort points in P around o, taking the direction of v as first angle
	sort(P.begin(), P.end(), [&](const point & a, const point & b){
		return point((a - o).half(v), 0) < point((b - o).half(v), (a - o).cross(b - o));
	});
}

struct coords{
    point3 o, dx, dy, dz;

    coords(point3 p, point3 q, point3 r) : o(p) {
        dx = (q - p).unit();
        dz = ((dx).cross(r - p)).unit();
        dy = dz.cross(dx);
    }

    coords(point3 p, point3 q, point3 r, point3 s) : 
     o(p), dx(q - p), dy(r - p), dz(s - p) {}
    //Return a point in 2D
    point pos2d(point3 p) {return point((p - o).dot(dx), (p - o).dot(dy));}
    point3 pos3d(point3 p) {return point3((p - o).dot(dx), (p - o).dot(dy), (p - o).dot(dz));}
};

int sgn(ld x){
    if(ge(x, 0)) return 1;
    if(le(x, 0)) return -1;
    return 0;
}

vector<pair<point, int>> convexHull(vector<pair<point, int>> P){
	sort(P.begin(), P.end());
	vector<pair<point, int>> L, U;
	for(int i = 0; i < P.size(); i++){
		while(L.size() >= 2 && leq((L[L.size() - 2].fi - P[i].fi).cross(L[L.size() - 1].fi - P[i].fi), 0)){
			L.pop_back();
		}
		L.push_back(P[i]);
	}
	for(int i = P.size() - 1; i >= 0; i--){
		while(U.size() >= 2 && leq((U[U.size() - 2].fi - P[i].fi).cross(U[U.size() - 1].fi - P[i].fi), 0)){
			U.pop_back();
		}
		U.push_back(P[i]);
	}
	L.pop_back();
	U.pop_back();
	L.insert(L.end(), U.begin(), U.end());
	return L;
}

bool pointInLine(const point & a, const point & v, const point & p){
	//line a+tv, point p
	return eq((p - a).cross(v), 0);
}

bool pointInSegment(const point & a, const point & b, const point & p){
	//segment ab, point p
	return pointInLine(a, b - a, p) && leq((a - p).dot(b - p), 0);
}

bool pointInPerimeter(const vector<point> & P, const point & p){
	int n = P.size();
	for(int i = 0; i < n; i++){
		if(pointInSegment(P[i], P[(i + 1) % n], p)){
			return true;
		}
	}
	return false;
}

bool crossesRay(const point & a, const point & b, const point & p){
	return (geq(b.y, p.y) - geq(a.y, p.y)) * sgn((a - p).cross(b - p)) > 0;
}

int pointInPolygon(vector<pair<point, int>> S, point p){
    vector<point> P;
    for(auto e : S) P.pb(e.fi);
	if(pointInPerimeter(P, p)){
		return 1; //point in the perimeter
	}
	int n = P.size();
	int rays = 0;
	for(int i = 0; i < n; i++){
		rays += crossesRay(P[i], P[(i + 1) % n], p);
	}
	return rays & 1; //0: point outside, 1: point inside
}

//Auxiliar points
point3 zero;

ld volume(vector<vector<point3>> P){
    int n = P.size();
    ld ans = 0;
    for(int i = 0; i < n; i++){
        point3 res;
        int m = P[i].size();
        for(int j = 0; j < m; j++){
            res += (P[i][j]).cross(P[i][(j + 1) % m]);
        }
        ans += (P[i][0]).dot(res);
    }
    return abs(ans) / 6.0;
}

point3 closestPointLs(const line & l1, const line & l2){
    point3 n = (l2.q).cross((l1.q).cross(l2.q));
    return l1.p + l1.q * (l2.p - l1.p).dot(n) / (l1.q).dot(n);
}

int infointersectSegments(point3 a, point3 b, point3 c, point3 d){
    line l1(a, b - a), l2(c, d - c);
    if((l1.q).cross(l2.q) == zero) return 0;
    point3 it = closestPointLs(l1, l2);
    if(it == a || it == b || it == c || it == d) return -1;
    if(it > min(a, b) && it < max(a, b) && it > min(c, d) && it < max(c, d)) return 1;
    return 0;
}

//Assuming that they intersect and not parallel, (pl.n).cross(ql.n) != {0,0,0}
line intersectionPlanes(const plane & pl, const plane & ql){
    point3 q = (pl.n).cross(ql.n);
    point3 p = (ql.n * pl.d - pl.n * ql.d).cross(q) / q.norm();
    return line(p, q);
}

bool point3InLine(const point3 & a, const point3 & v, const point3 & p){
	return eq(line(a, v).dist(p), 0);
}

bool point3InSegment(const point3 & a, const point3 & b, const point3 & p){
    if(p == a || p == b) return 1;
	return point3InLine(a, b - a, p) && p >= min(a, b) && p <= max(a, b);
}

int intersectSegmentPlaneInfo(const point3 & a, const point3 & b, const plane & p){
    //if(leq(p.dist(a), 0) && leq(p.dist(b), 0)) return -1;
    point3 d = b - a;
    if(eq(d.dot(p.n), 0)){
        if(eq(p.dist(a), 0)) return -1;
        else return 0;
    }else{
        if(point3InSegment(a, b, line(a, b - a).inter(p))) return 1;
    }
    return 0;
}

vector<vector<point3>> reorient(vector<vector<point3>> P){
    int n = P.size();
    vector<vector<pair<int, bool>>> g(n);
    map<pair<point3, point3>, int> es;
    for(int i = 0; i < n; i++){
        int m = P[i].size();
        for(int j = 0; j < m; j++){
            point3 a = P[i][j], b = P[i][(j + 1) % m];
            if(!es.size()){
                es[{a, b}] = i;
            }else if(es.count({a, b})){
                int v = es[{a, b}];
                g[i].pb({v, true});
                g[v].pb({i, true});
            }else if(es.count({b, a})){
                int v = es[{b, a}];
                g[i].pb({v, false});
                g[v].pb({i, false});
            }else es[{a, b}] = i;
        }
    }
    vector<bool> vis(n, false), flip(n, false);
    flip[0] = false;
    queue<int> q;
    q.push(0);
    while(q.size()){
        int u = q.front();
        q.pop();
        for(auto e : g[u]){
            if(!vis[e.fi]){
                vis[e.fi] = true;
                flip[e.fi] = (flip[u] ^ e.se);
                q.push(e.fi);
            }
        }
    }
    for(int i = 0; i < n; i++){
        if(flip[i]) reverse(all(P[i]));
    }
    return P;
}

int vls(bool a, bool b){
    if(a && b) return 0;
    if(!a && !b) return 2;
    if(a) return 1;
    return -1;
}

int sidePolyhedroPlane(vector<vector<point3>> P, plane p){
    int n = P.size(); bool sup = false, infr = false;
    for(int i = 0; i < n; i++){
        for(int l = 0; l < P[i].size(); l++){
            auto u = p.side(P[i][l]);
            if(ge(u, 0)) sup = true;
            if(le(u, 0)) infr = true;
        }
    }
    return vls(sup, infr);
}

int sideFacePlane(vector<point3> F, plane p){
    int n = F.size(); bool sup = false, infr = false;
    for(int i = 0; i < n; i++){
        auto u = p.side(F[i]);
        if(ge(u, 0)) sup = true;
        if(le(u, 0)) infr = true;
    }
    return vls(sup, infr);
}

vector<pair<point, int>> createFaces3(vector<vector<int>> L, vector<pair<point, int>> vp, vector<int> & flags, int i0, int i1){
    //BFS
    queue<pair<point, int>> q;
    stack<pair<point, int>> s;
    s.push(vp[i0]); q.push(vp[i0]);
    while(q.size()){
        int v = (q.front()).se; q.pop();
        if(v == -1) return {};
        if(flags[v]) break;
        flags[v] = true;
        int i = i1, next = -1; ld an = 0, alp = (vp[i].fi - vp[v].fi).ang();
        for(int j = 0; j < L[v].size(); j++){
            if(L[v][j] == v || L[v][j] == i) continue;
            ld se = (vp[L[v][j]].fi - vp[v].fi).ang(); se -= alp;
            if(le(se, 0)) se += pi * 2;
            if(ge(se, an)){
                an = se;
                next = L[v][j];
            }
        }
        q.push(vp[next]); s.push(vp[next]); i1 = v;
    }
    if(!s.size()) return {};
    auto v = s.top(); s.pop();
    vector<pair<point, int>> res;
    res.pb(v);
    while(s.size()){
        auto u = s.top();
        if(u.se == v.se) break;
        res.pb(u); s.pop();
    }
    return res;
}

bool polygonInPolygon(vector<pair<point, int>> A, vector<pair<point, int>> B){
    for(int i = 0; i < A.size(); i++){
        if(!pointInPolygon(B, A[i].fi)) return false;
    }
    return true;
}

vector<pair<point, int>> createFaces2(vector<vector<int>> L, vector<pair<point, int>> vp, vector<int> & flags, int i){
    auto f = flags, g = flags;
    auto ft = createFaces3(L, vp, f, i, L[i][0]); 
    auto tf = createFaces3(L, vp, g, L[i][0], i);
    if(ft.size() < 3) return {};
    if(tf.size() < 3) return {};
    if(convexHull(ft).size() < 3) return {};
    if(convexHull(tf).size() < 3) return {};
    if(polygonInPolygon(tf, ft)){
        for(int i = 0; i < flags.size(); i++) flags[i] = f[i];
        return ft;   
    }else if(polygonInPolygon(ft, tf)){
    	for(int i = 0; i < flags.size(); i++) flags[i] = g[i];
    	return tf;
    }
    return {};
}

vector<vector<pair<point, int>>> createFaces1(vector<vector<int>> L, vector<pair<point, int>> vp){
    vector<vector<pair<point, int>>> ans;
    int n = vp.size();
    vector<int> flags(n);
    for(int i = 0; i < n; i++){
        if(!flags[i] && L[i].size() != 0){
            auto v = createFaces2(L, vp, flags, i);
            if(v.size() > 2) ans.pb(v);
        }
    }
    return ans;
}

vector<vector<point3>> createFaces(vector<vector<int>> L, vector<point3> pts){
     vector<vector<point3>> ans;
    int n = pts.size();
    if(n < 3) return {};
    vector<bool> flags(n);
    int j = 2;
    while((pts[1] - pts[0]).cross(pts[j] - pts[0]) == zero && j < n) j++;
    if(j == n) return {};
    coords c(pts[0], pts[1], pts[j]);
    vector<pair<point, int>> vp;
    for(int i = 0; i < n; i++) vp.pb({c.pos2d(pts[i]), i});
    auto ft = createFaces1(L, vp);
    for(int i = 0; i < ft.size(); i++){
        vector<point3> v;
        for(auto e : ft[i]) v.pb(pts[e.se]);
        ans.pb(v);
    }
    return ans;
}

bool connectFaces(vector<point3> A, vector<point3> B){
    int n = A.size(), m = B.size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            point3 a = min(A[i], A[(i + 1) % n]), b = max(A[i], A[(i + 1) % n]), c = min(B[j], B[(j + 1) % m]), d = max(B[j], B[(j + 1) % m]);
            if(a == c && b == d) return true;
            if(!((a - b).cross(c - d) == zero)) continue;
            if(neq(line(a, b - a).dist(c), 0)) continue;
            if(a >= d || b <= c) continue;
            return true;
        }
    }
    return false;
}

void createPolyhedrosDFS(vector<vector<int>> & L, vector<bool> & f, int in, vector<int> & res){
    res.pb(in); f[in] = true;
    for(int i = 0; i < L[in].size(); i++){
        if(!f[L[in][i]]){
            createPolyhedrosDFS(L, f, L[in][i], res);
        }
    }
    return;
}

vector<vector<vector<point3>>> createPolyhedros(vector<vector<point3>> F){
    int n = F.size();
    vector<vector<vector<point3>>> Ps;
    if(n < 3) return {};
    vector<vector<int>> L(n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i == j) continue;
            if(connectFaces(F[i], F[j])) L[i].pb(j);
        }
    }
    vector<bool> flags(n, false);
    vector<vector<int>> sp;
    for(int i = 0; i < n; i++){
        if(!flags[i]){
            vector<int> ts;
            createPolyhedrosDFS(L, flags, i, ts);
            if((int) ts.size() > 3) sp.pb(ts);
        }
    }
    for(int i = 0; i < sp.size(); i++){
        vector<vector<point3>> ts;
        for(int j = 0; j < sp[i].size(); j++){
            ts.pb(F[sp[i][j]]);
        }
        if(ts.size() > 3) Ps.pb(ts);
    }
    return Ps;
}

vector<point3> clearFace(vector<point3> P, plane p){
    vector<point3> ans;
    int n = P.size();
    vector<bool> fl(n);
    if(P[0] < P[n - 1]){
        if(eq(p.side(P[0]), 0) && eq(p.side(P[1]), 0)) if(P[1] > P[0]) fl[0] = true;
        if(eq(p.side(P[n - 1]), 0) && eq(p.side(P[n - 2]), 0)) if(P[n - 2] < P[n - 1]) fl[n - 1] = true;
    }else{
        if(eq(p.side(P[0]), 0) && eq(p.side(P[1]), 0)) if(P[1] < P[0]) fl[0] = true;
        if(eq(p.side(P[n - 1]), 0) && eq(p.side(P[n - 2]), 0)) if(P[n - 2] > P[n - 1]) fl[n - 1] = true;
    }
    for(int i = 0; i < n; i++){
        if(!fl[i]) ans.pb(P[i]);
    }
    return ans;
}

bool verify1(vector<point3> A){
    int n = A.size(), j = 2;
    while((A[j] - A[0]).cross(A[1] - A[0]) == zero && j < n) j++;
    if(j == n) return false;
    return true;
}

vector<vector<point3>> cutFace(vector<point3> F, plane p){
    int n = F.size();
    if(n < 3) return {};
    vector<vector<point3>> ans;
    int i = 0;
    while(geq(p.side(F[i]), 0) && i < n) i++;
    if(i == n) return {};
    int fin = i;
    vector<point3> res;
    bool flag = false;
    do{
        if(!flag){
            if(intersectSegmentPlaneInfo(F[i], F[(i + 1) % n], p) == 1){
                point3 pt = line(F[i], F[(i + 1) % n] - F[i]).inter(p);
                if(pt != F[(i + 1) % n]){
                    res.pb(pt);
                    flag = true;
                }
            }else if(intersectSegmentPlaneInfo(F[i], F[(i + 1) % n], p) == -1){
                res.pb(F[i]);
                flag = true;
            }
        }else{
            if(le(p.side(F[i]), 0)){
                if(res.size() > 2) ans.pb(res);
                res.clear();
                flag = false;
            }
            if(geq(p.side(F[i]), 0)) res.pb(F[i]);
            if(intersectSegmentPlaneInfo(F[i], F[(i + 1) % n], p) == 1){
                auto pt = line(F[i], F[(i + 1) % n] - F[i]).inter(p);
                if(pt != F[i] && pt != F[(i + 1) % n]){
                    res.pb(pt);
                    if(res.size() > 2) ans.pb(res);
                    res.clear();
                    flag = false;
                }
            }
        }
        i = (i + 1) % n;
    }while(i != fin);
    if((int)res.size() > 2) ans.pb(res);
    vector<vector<point3>> ans1;
    for(int i = 0; i < ans.size(); i++){
        auto r = clearFace(ans[i], p);
        if(r.size() > 2) if(verify1(r)) ans1.pb(r);
    }
    return ans1;
}

vector<vector<vector<point3>>> cutPolyhedro(vector<vector<point3>> P, plane p){
    int n = P.size();
    if(n < 4) return {};
    vector<vector<point3>> ans, auxln;
    if(sidePolyhedroPlane(P, p) == 1){
        vector<vector<vector<point3>>> ans1; ans1.pb(P); return ans1;    
    }
    for(int i = 0; i < n; i++){
        if(P[i].size() < 2) continue;
        auto ri = sideFacePlane(P[i], p);
        if(ri != 0){
            if(ri == 2){
                vector<point3> ax;
                ax.pb(P[i][P[i].size() - 1]);
                for(auto e : P[i]) ax.pb(e);
                auxln.pb(ax);
                continue;
            }else{
                int j = 0, m = P[i].size(); 
                while(eq(p.side(P[i][j]), 0) && j < m) j++;
                if(j == m) j = 0;
                vector<point3> ax;
                int fn = j; bool fl = false;
                do{
                    if(!fl){
                        if(eq(p.side(P[i][j]), 0)){
                            ax.pb(P[i][j]); fl = true;
                        }
                    }else{
                        if(eq(p.side(P[i][j]), 0)){
                            ax.pb(P[i][j]);
                        }else{
                            if(ax.size() > 1) auxln.pb(ax);
                            ax.clear(); fl = false;
                        }
                    }
                j = (j + 1) % m;
                }while(j != fn);
                if(ax.size() > 1) auxln.pb(ax);
                if(ri == 1) ans.pb(P[i]);
                continue;
            }
        }
        vector<vector<point3>> res = cutFace(P[i], p);
        for(int j = 0; j < res.size(); j++){
            ans.pb(res[j]);
        }
        for(int l = 0; l < res.size(); l++){
            int j = 0, m = res[l].size(); 
            while(eq(p.side(res[l][j]), 0) && j < m) j++;
            if(j == m) j = 0;
            vector<point3> ax;
            int fn = j; bool fl = false;
            do{
                if(!fl){
                    if(eq(p.side(res[l][j]), 0)){
                        ax.pb(res[l][j]); fl = true;
                    }
                }else{
                    if(eq(p.side(res[l][j]), 0)){
                        ax.pb(res[l][j]);
                    }else{
                        if(ax.size() > 1) auxln.pb(ax);
                        ax.clear(); fl = false;
                    }
                }
            j = (j + 1) % m;
            }while(j != fn);
            if(ax.size() > 1) auxln.pb(ax);
        }
    }
    vector<point3> pts;
    for(auto e : auxln) for(auto f : e) pts.pb(f);
    if(pts.size()){
        sort(all(pts));
        pts.erase(unique(all(pts)), pts.end());
        if(pts.size() < 3) return createPolyhedros(ans);
        vector<vector<int>> L;
        L.resize((int) pts.size());
        for(int i = 0; i < auxln.size(); i++){
            for(int j = 0; j < auxln[i].size() - 1; j++){
                int it = lower_bound(all(pts), auxln[i][j]) - pts.begin(), at = lower_bound(all(pts), auxln[i][j + 1]) - pts.begin();
                L[it].pb(at); L[at].pb(it);
            }
        }
        for(int i = 0; i < pts.size(); i++){
            if(!L[i].size()) continue;
            sort(all(L[i]));
            L[i].erase(unique(all(L[i])), L[i].end());
        }
        auto sup = createFaces(L, pts);
        for(auto e : sup) ans.pb(e);
    }
    return createPolyhedros(ans);
}

vector<vector<point3>> inicialize(){
    vector<vector<point3>> faces;
    faces.pb({point3(-Inf, -Inf, -Inf), point3(Inf, -Inf, -Inf), point3(Inf, Inf, -Inf), point3(-Inf, Inf, -Inf)});
    faces.pb({point3(-Inf, -Inf, -Inf), point3(-Inf, -Inf, Inf), point3(-Inf, Inf, Inf), point3(-Inf, Inf, -Inf)});
    faces.pb({point3(-Inf, -Inf, -Inf), point3(Inf, -Inf, -Inf), point3(Inf, -Inf, Inf), point3(-Inf, -Inf, Inf)});
    faces.pb({point3(-Inf, Inf, -Inf), point3(Inf, Inf, -Inf), point3(Inf, Inf, Inf), point3(-Inf, Inf, Inf)});
    faces.pb({point3(Inf, -Inf, -Inf), point3(Inf, Inf, -Inf), point3(Inf, Inf, Inf), point3(Inf, -Inf, Inf)});
    faces.pb({point3(-Inf, -Inf, Inf), point3(Inf, -Inf, Inf), point3(Inf, Inf, Inf), point3(-Inf, Inf, Inf)});
    return faces;
}

int main(){
    //ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int n; cin>>n;
    vector<vector<point3>> A;
    for(int i = 0; i < n; i++){
        int k; cin>>k;
        vector<point3> C(k);
        for(int j = 0; j < k; j++){
            int a, b, c; cin>>a>>b>>c;
            C[j] = point3(a, b, c);
        }
        A.pb(C);
    }
    A = reorient(A);
    if(ge(volume(A), 0)) for(auto & e : A) reverse(all(e));
    point3 p0; ld d;
    cin>>p0.x>>p0.y>>p0.z>>d;
    plane pl(p0, d);
    auto B = cutPolyhedro(A, pl);
    cout<<B.size()<<'\n';
    for(int i = 0; i < B.size(); i++){
        cout<<B[i].size();
        for(int j = 0; j < B[i].size(); j++){
            cout<<B[i][j].size();
            for(int l = 0; l < B[i][j].size(); l++){
                cout<<B[i][j][l]<<" ";
            }
            cout<<'\n';
        }
        cout<<'\n';
    }
}
