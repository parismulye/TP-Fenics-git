// Gmsh project created on Mon Sep 14 09:33:57 2015
SetFactory("OpenCASCADE");
a = 1.;
b = 1.;
R = 0.2;

//seed size:
d = 0.03;

Point(1) = {0, 0, 0, d};
Point(2) = {a, 0, 0, d};
Point(3) = {a, b, 0, d};
Point(4) = {0, b, 0, d};
Point(5) = {R,0,0,d};
Point(6) = {a-R,0,0,d};
Point(7) = {a,R,0,d};
Point(8) = {a,b-R,0,d};
Point(9) = {a-R,b,0,d};
Point(10) = {R,b,0,d};
Point(11) = {0,b-R,0,d};
Point(12) = {0,R,0,d};


Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 2};
Line(4) = {2, 7};
Line(5) = {7, 8};
Line(6) = {8, 3};
Line(7) = {3, 9};
Line(8) = {9, 10};
Line(9) = {10, 4};
Line(10) = {4, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};
Circle(13) = {5, 1, 12};
Circle(14) = {7, 2, 6};
Circle(15) = {9, 3, 8};
Circle(16) = {11, 4, 10};
Ellipse(17) = {0.5, 0.5, 0, 0.4, 0.2, 0, 2*Pi};



Line Loop(17) = {11, -13, 2, -14, 5, -15, 8, -16};
// added circle in center
Line Loop(27) = {17};
Plane Surface(18) = {17, 27};
Line Loop(19) = {1, 13, 12};
Plane Surface(20) = {19};
Line Loop(21) = {3, 4, 14};
Plane Surface(22) = {21};
Line Loop(23) = {7, 15, 6};
Plane Surface(24) = {23};
Line Loop(25) = {10, 16, 9};
Plane Surface(26) = {25};

// added circle in center surface
Line Loop(30) = {17};
Plane Surface(27) = {30};
Transfinite Line{17} = 3*Pi*R/d;


Physical Line(1) = {10,11,12};
Physical Line(2) = {1,2,3};
Physical Line(3) = {4,5,6};
Physical Line(4) = {7,8,9};

Physical Surface(0) = {18};               // main domain has ID=0
Physical Surface(1) = {20,22,24,26,27};   // central inclusion has ID=1

