#include <bits/stdc++.h>
#define T double
#define pi acos(-1)
#define mk make_pair
#define pb push_back

using namespace std;

struct pt{
	T x, y;	//define T as int, ll, or double depending on problem
	pt operator+(pt p){return {x + p.x, y + p.y};}
	pt operator-(pt p){return {x - p.x, y - p.y};}
	pt operator*(T d){return {x * d, y * d};}
	pt operator/(T d){return {x / d, y / d};} // only for doubles

};

bool operator==(pt a, pt b){return a.x == b.x and a.y == b.y;}
bool operator!=(pt a, pt b){return !(a == b);}
bool operator<(pt a, pt b){return mk(a.y, a.x) < mk(b.y, b.x);}	//for use in sets
T sq(pt p) {return p.x*p.x + p.y*p.y;}	//square
double abs(pt p) {return sqrt(sq(p));}	//magnitude
T dot(pt v, pt w){return v.x * w.x + v.y * w.y;}
T cross(pt v, pt w){return v.x * w.y - v.y * w.x;}
bool isPerp(pt v, pt w) {return dot(v,w) == 0;}
T orient(pt a, pt b, pt c) {return cross(b-a,c-a);} //going from a to b to c; positive if left, negative if right, 0 if collinear
double angle(pt v, pt w){	//angle between 2 vectors. [0, pi]
	double cosTheta = dot(v, w) / abs(v) / abs(w);
	return acos(max(-1.0, min(1.0, cosTheta)));
}
pt rot(pt p, double a){	return {p.x*cos(a) - p.y*sin(a), p.x*sin(a) + p.y*cos(a)};}
T area2(pt a, pt b, pt c) { return cross(a,b) + cross(b,c) + cross(c,a); } //triangle area
double areaPolygon(vector<pt> p){
	double area = 0.0;
	int s = p.size();
	for(int i = 0; i < s; i++)
		area += cross(p[i], p[(i + 1) % s]);
	return abs(area) / 2.0;
}
ostream& operator<<(ostream& os, pt p) {	//for debugging
	return os << "("<< p.x << "," << p.y << ")";
}

//LINES
struct line{
	pt v; T c;	
	line(pt v, T c) : v(v), c(c) {} //from direction v, offset c
	line(T a, T b, T c) : v({b, -a}), c(c) {} //from ax+by=c
	line(pt p, pt q) : v(q - p), c(cross(v, p)) {} //from ax+by=c

	T side(pt p) {return cross(v, p) - c;}
	double dist(pt p) {return abs(side(p)) / abs(v);}
	double sqDist(pt p) {return side(p)*side(p) / (double)sq(v);}
	line translate(pt t) {return {v, c + cross(v,t)};}
	bool cmpProj(pt p, pt q) {
		return dot(v,p) < dot(v,q);
	}

};

bool inter(line l1, line l2, pt &out) {	//do l1 and l2 interesect? 
	T d = cross(l1.v, l2.v);
	if (d == 0) return false;
	out = (l2.v * l1.c - l1.v * l2.c) / d; // requires floating-point coordinates
	return true;
}
//SEGMENTS
bool inDisk(pt a, pt b, pt p) { //is p in the circle with diameter ab
	return dot(a-p, b-p) <= 0;
}
bool onSegment(pt a, pt b, pt p) {	//is p on ab?
	return orient(a,b,p) == 0 && inDisk(a,b,p);
}
bool properInter(pt a, pt b, pt c, pt d, pt &out) { //determines if ab and cd interesect in the middle
	double oa = orient(c,d,a),
	ob = orient(c,d,b),
	oc = orient(a,b,c),
	od = orient(a,b,d);
	// Proper intersection exists iff opposite signs
	if (oa*ob < 0 && oc*od < 0) {
		out = (a*ob - b*oa) / (ob-oa);
		return true;
	}
	return false;
}
set<pt> inters(pt a, pt b, pt c, pt d) {	//returns 0, 1, or 2 points describing empty, single point, or segment
	pt out;
	if (properInter(a,b,c,d,out)) return {out};
	set<pt> s;
	if (onSegment(c,d,a)) s.insert(a);
	if (onSegment(c,d,b)) s.insert(b);
	if (onSegment(a,b,c)) s.insert(c);
	if (onSegment(a,b,d)) s.insert(d);
	return s;
}

double segPoint(pt a, pt b, pt p) { //dist between segment AB and p
	if (a != b) {
		line l(a,b);
		if (l.cmpProj(a,p) && l.cmpProj(p,b)) // if closest to projection
			return l.dist(p); // output distance to line
	}
	return min(abs(p-a), abs(p-b)); // otherwise distance to A or B
}

double segSeg(pt a, pt b, pt c, pt d) { //dist between segments
	pt dummy;
	if (properInter(a,b,c,d,dummy))
		return 0;
	return min({segPoint(a,b,c), segPoint(a,b,d),
		segPoint(c,d,a), segPoint(c,d,b)});
}

//Check if a point is in a polygon
//should be counterclockwise order?

// true if P at least as high as A (blue part)
bool above(pt a, pt p) {
	return p.y >= a.y;
}
// check if [PQ] crosses ray from A
bool crossesRay(pt a, pt p, pt q) {
	return (above(a,q) - above(a,p)) * orient(a,p,q) > 0;
}
// if strict, returns false when A is on the boundary
bool inPolygon(vector<pt> p, pt a, bool strict = true) {
	int numCrossings = 0;
	for (int i = 0, n = p.size(); i < n; i++) {
		if (onSegment(p[i], p[(i+1)%n], a))
			return !strict;
		numCrossings += crossesRay(a, p[i], p[(i+1)%n]);
	}
	return numCrossings & 1; // inside if odd number of crossings
}

vector<pt> convexhull(vector<pt> p){	//counterclockwise, no collinear points
	sort(p.begin(), p.end());
	p.erase(unique(p.begin(), p.end()), p.end());
	vector<pt> up, dn;
	for(pt i : p){
		while(up.size() > 1 and orient(up[up.size() - 2], up.back(), i) >= 0) up.pop_back();
		while(dn.size() > 1 and orient(dn[dn.size() - 2], dn.back(), i) <= 0) dn.pop_back();
		up.pb(i);
		dn.pb(i);
	}
	for(int i = (int) up.size() - 2; i >= 1; i--)
		dn.pb(up[i]);
	return dn;
}

int main(){
	pt a{3, 4}, b{2, -1};
	cout << a / 2 << " " << a - b << endl;
}