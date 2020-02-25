// Gmsh project created on Mon Sep 14 09:33:57 2015
SetFactory("OpenCASCADE");
L = 1.;
H = 1.;
R = 0.3;

//seed size:
d = 0.03;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};
Point(3) = {L, H, 0, d};
Point(4) = {0, H, 0, d};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Circle(5) = {L/2, H/2, 0, R, 0, 2*Pi};


Line Loop(1) = {1,2,3,4};
//added circle in center
Line Loop(2) = {5};
Plane Surface(1) = {1, 2};
Plane Surface(2) = {2};


Transfinite Line{5} = (2*Pi*R)/d;

Physical Line(1) = {4};
Physical Line(2) = {1};
Physical Line(3) = {2};
Physical Line(4) = {3};

Physical Surface(0) = {1};
Physical Surface(1) = {2};
