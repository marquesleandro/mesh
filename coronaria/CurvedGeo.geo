// Gmsh project created on Wed Dec 21 23:16:33 2016

lc = DefineNumber[ 0.02, Name "Parameters/lc" ]; 

// removing 1.0 from y coordinate for symmetry, since 
// normalV should be at Y=0.0 and everything else above. 
// This number makes the average channel diameter equal to 1.
xMax = 10.0000;
symY = 1.0;
Point(1)={0.000,symY-0,0,lc}; 
Point(2)={1.725,symY-0,0,lc}; 
Point(3)={8.275,symY-0,0,lc}; 
Point(4)={xMax,symY-0,0,lc}; 
Point(5)={xMax,symY-1,0,lc}; 
Point(6)={0,symY-1,0,lc}; 
Point(7)={5,symY-0.6,0,lc};
Point(8)={3,symY-0,0,lc}; 
Point(9)={7,symY-0,0,lc}; 

BSpline(1) = {2, 8, 7, 9, 3};
Line(2) = {1, 2};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line Loop(11) = {2, 1, 3, 4, 5, 6};
Plane Surface(11) = {11};

Physical Line("dirichlet1 noslip") = {2,1,3};
Physical Line("neumann1 outflow") = {4};
Physical Line("dirichlet1 inflow") = {6};
Physical Line("neumann1 symmetry") = {5};
Physical Line("dirichlet2 noslip") = {2,1,3};
Physical Line("neumann2 outflow") = {4};
Physical Line("dirichlet2 inflow") = {6};
Physical Line("dirichlet2 symmetry") = {5};
Physical Line("dirichlet3 bottom") = {5};
Physical Line("neumann3 outflow") = {4};
Physical Line("neumann3 inflow") = {6};
Physical Line("dirichlet3 top") = {2,1,3};
Physical Line("neumann4 bottom") = {5};
Physical Line("neumann4 outflow") = {3,4};
Physical Line("neumann4 inflow") = {2,6};
Physical Line("dirichlet4 top") = {1};

Physical Surface(11) = {11};

