// Gmsh project created on Sun Feb 23 19:44:01 2020
SetFactory("OpenCASCADE");
//+
r = 0.1330;
//+
h = 0.0337;
d = 0.0505;
Circle(1) = {d+r-h, d+r-h, 0, r, 0, 2*Pi};
//+
Circle(2) = {d+r+r+d+r, d+r-h, 0, r, 0, 2*Pi};
//+
Circle(3) = {d+r+r+d+r+r+d+r+h, d+r-h, 0, r, 0, 2*Pi};
//+
Circle(4) = {d+r-h, d+r+r+d+r, 0, r, 0, 2*Pi};
//+
Circle(5) = {d+r+r+d+r, d+r+r+d+r, 0, r, 0, 2*Pi};
//+
Circle(6) = {d+r+r+d+r+r+d+r+h, d+r+r+d+r, 0, r, 0, 2*Pi};

Circle(7) = {d+r-h, d+r+r+d+r+r+d+r+h, 0, r, 0, 2*Pi};
//+
Circle(8) = {d+r+r+d+r, d+r+r+d+r+r+d+r+h, 0, r, 0, 2*Pi};
//+
Circle(9) = {d+r+r+d+r+r+d+r+h, d+r+r+d+r+r+d+r+h, 0, r, 0, 2*Pi};

//+
Point(13) = {0, 0, 0, 0.02};
//+
Point(14) = {1, 0, 0, 0.02};
//+
Point(15) = {1, 1, 0, 0.02};
//+
Point(16) = {0, 1, 0, 0.02};
//+
Line(17) = {16, 13};
//+
Line(18) = {13, 14};
//+
Line(19) = {14, 15};
//+
Line(20) = {15, 16};
//+
Line Loop(1) = {17, 18, 19, 20};
//+
//+
Line Loop(2) = {17, 18, 19, 20};
//+
Line Loop(3) = {7};
//+
Line Loop(4) = {8};
//+
Line Loop(5) = {9};
//+
Line Loop(6) = {4};
//+
Line Loop(7) = {5};
//+
Line Loop(8) = {6};
//+
Line Loop(9) = {1};
//+
Line Loop(10) = {2};
//+
Line Loop(11) = {3};
//+
Plane Surface(1) = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
//+
Line Loop(12) = {7};
//+
Plane Surface(2) = {12};
//+
Line Loop(13) = {8};
//+
Plane Surface(3) = {13};
//+
Line Loop(14) = {9};
//+
Plane Surface(4) = {14};
//+
Line Loop(15) = {4};
//+
Plane Surface(5) = {15};
//+
Line Loop(16) = {5};
//+
Plane Surface(6) = {16};
//+
Line Loop(17) = {6};
//+
Plane Surface(7) = {17};
//+
Line Loop(18) = {2};
//+
Plane Surface(8) = {18};
//+
Line Loop(19) = {1};
//+
Plane Surface(9) = {19};
//+
Line Loop(20) = {3};
//+
Plane Surface(10) = {20};
//+
Physical Line(1) = {17};
//+
Physical Line(2) = {18};
//+
Physical Line(3) = {19};
//+
Physical Line(4) = {20};
//+
Physical Surface(0) = {1};
//+
Physical Surface(1) = {2, 3, 4, 5, 6, 7, 9, 8, 10};
//+
Transfinite Line {7, 8, 9, 4, 5, 6, 1, 2, 3} = 60 Using Progression 1;
