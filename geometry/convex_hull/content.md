```sh
struct point {
    ld x, y;
    point(ld xloc, ld yloc) : x(xloc), y(yloc) {}
    point() {}
    point& operator= (const point& other) {
        x = other.x, y = other.y;
        return *this;
    }
    bool operator== (const point& other) const {
        return abs(other.x - x) < prec && abs(other.y - y) < prec;
    }
    bool operator!= (const point& other) const {
        return !(abs(other.x - x) < prec && abs(other.y - y) < prec);
    }
    bool operator< (const point& other) const {
        return (x < other.x ? true : (x == other.x && y < other.y));
    }
};
```
```sh
struct vect {
    ld i, j;
};
```
```sh
struct segment {
    point p1, p2;
    segment(point a, point b) : p1(a), p2(b) {}
    segment() {}
};
```
```sh
ld crossproduct(point a, point b, point c) {
    vect ab, ac;
    ab.i = b.x - a.x;
    ab.j = b.y - a.y;
    ac.i = c.x - a.x;
    ac.j = c.y - a.y;
    return (ab.i * ac.j - ab.j * ac.i);
}
```
```sh
ld dotproduct(point a, point b, point c) {
    vect ab, bc;
    ab.i = b.x - a.x;
    ab.j = b.y - a.y;
    bc.i = c.x - b.x;
    bc.j = c.y - b.y;
    return (ab.i * bc.i + ab.j * bc.i);
}
```
```sh
int orientation(point p, point q, point r) {
    int val = (int) crossproduct(p, q, r);
    if(val == 0) return COLINEAR;
    return (val > 0) ? CW : CCW;
}
```
```sh
bool onsegment(point p, segment s) {
    return (p.x <= max(s.p1.x, s.p2.x) && p.x >= min(s.p1.x, s.p2.x) &&
            p.y <= max(s.p1.y, s.p2.y) && p.y >= min(s.p1.y, s.p2.y));
}
```
```sh
ld pointdist(point a, point b) {
    return sqrt(pow(a.x-b.x, 2) + pow(a.y-b.y, 2));
}
```
```sh
ld pointsquareddist(point a, point b) {
    return pow(a.x-b.x, 2) + pow(a.y-b.y, 2);
}
```
```sh
vector<point> intersect(segment s1, segment s2) {
    vector<point> res;
    point a = s1.p1, b = s1.p2, c = s2.p1, d = s2.p2;
    if(orientation(a,b,c) == 0 && orientation(a,b,d) == 0 &&
       orientation(c,d,a) == 0 && orientation(c,d,b) == 0) {
        point min_s1 = min(a,b), max_s1 = max(a,b);
        point min_s2 = min(c,d), max_s2 = max(c,d);
        if(min_s1 < min_s2) {
            if(max_s1 < min_s2) return res;
        }
        else if(min_s2 < min_s1 && max_s2 < min_s1) return res;
        point start = max(min_s1,min_s2), end = min(max_s1,max_s2);
        if(start == end) res.pb(start);
        else {
            res.pb(min(start,end));
            res.pb(max(start,end));
        }
        return res;
    }
    ld x1 = (b.x-a.x);
    ld y1 = (b.y-a.y);
    ld x2 = (d.x-c.x);
    ld y2 = (d.y-c.y);
    ld u1 = (-y1 * (a.x-c.x) + x1 * (a.y-c.y)) / (-x2 * y1 + x1 * y2);
    ld u2 = (x2 * (a.y-c.y) - y2 * (a.x-c.x)) / (-x2 * y1 + x1 * y2);
    if(u1 >= 0 && u1 <= 1 && u2 >= 0 && u2 <= 1) {
        res.pb(point((a.x + u2 * x1), (a.y + u2 * y1)));
    }
    return res;
}
```
```sh
ld linepointdist(segment s, point p, bool issegment) {
    if(s.p1 == s.p2) {
        if(p == s.p1) return 0;
        return pointdist(p, s.p1);
    }
    if(issegment) {
        if(dotproduct(s.p1, s.p2, p) > 0) return pointdist(s.p2, p);
        if(dotproduct(s.p2, s.p1, p) > 0) return pointdist(s.p1, p);
    }
    return abs(crossproduct(s.p1,s.p2,p) / pointdist(s.p1,s.p2));
}
```
```sh
vector<point> convexhull(vector<point>& X, bool onedge = false) {
    vector<point> hull;
    ll N = X.size();
    vector<bool> used(N, false);
    int p = 0;
    for(int i = 1; i < N; i++) if(X[i].x < X[p].x) p = i;
    int start = p;
    do {
        int n = -1;
        ld dist = onedge ? iinf : 0;
        for(int i = 0; i < N; i++) {
            if(i == p || used[i]) continue;
            if(n == -1) n = i;
            ld cross = crossproduct(X[p], X[i], X[n]);
            ld d = pointsquareddist(X[i], X[p]);
            if(cross < 0 || (cross == 0 &&
                ((onedge && d < dist) || (!onedge && d > dist)))) {
                n = i;
                dist = d;
            }
        }
        p = n;
        used[p] = true;
        hull.pb(X[p]);
    } while(start != p);
    return hull;
}
```
```sh
ld polyarea(vector<point>& p) {
    ld result = 0;
    for(int i = 0, j = 1, n = p.size(); i < n; i++, j = (i+1)%n)
        result += (p[i].x * p[j].y) - (p[i].y * p[j].x);
    return result / 2.0;
}
```
```sh
int insidepoly(vector<point> poly, point p) {
    bool inside = false;
    point outside(iinf, p.y);
    vector<point> intersection;
    for(uint i = 0, j = poly.size()-1; i < poly.size(); i++, j=i-1) {
        if(p == poly[i] || p == poly[j]) return ONEDGE;
        if(orientation(p,poly[i],poly[j]) == COLINEAR &&
           onsegment(p,segment(poly[i],poly[j]))) return ONEDGE;
        intersection = intersect(segment(p,outside), segment(poly[i],poly[j]));
        if(intersection.size() == 1) {
            if(poly[i] == intersection[0] && poly[j].y <= p.y) continue;
            if(poly[j] == intersection[0] && poly[i].y <= p.y) continue;
            inside = !inside;
        }
    }
    return inside ? INSIDE : OUTSIDE;
}
```