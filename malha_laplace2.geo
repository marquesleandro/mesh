cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 1, 0, 1};
Point(3) = {1, 0, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 3};
Line(2) = {3, 2};
Line(3) = {2, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Physical Line("sup1") = {1};
Physical Line("sup2") = {2};
Physical Line("sup3") = {3};
Physical Line("sup4") = {4};
Physical Surface(11) = {6};