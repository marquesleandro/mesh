lc = DefineNumber[ 0.02, Name "Parameters/lc" ]; 
// removing 2.0-0.605 from y coordinate for symmetry, since 
// normalV should be at Y=0.0 and everything else above. 
// This number makes the average channel diameter equal to 1.
xMax = 10.0000;
xFactor = 1.0;
symY = 2.0-0.605;
Point(1)={0.000000*xFactor,symY-0.273000,0,lc}; 
Point(2)={0.102400*xFactor,symY-0.273000,0,lc}; 
Point(3)={0.136500*xFactor,symY-0.238900,0,lc}; 
Point(4)={0.170600*xFactor,symY-0.204800,0,lc}; 
Point(5)={0.273000*xFactor,symY-0.204800,0,lc}; 
Point(6)={0.307200*xFactor,symY-0.170600,0,lc}; 
Point(7)={0.341300*xFactor,symY-0.102400,0,lc}; 
Point(8)={0.375400*xFactor,symY-0.102400,0,lc}; 
Point(9)={0.409600*xFactor,symY-0.068300,0,lc}; 
Point(10)={0.580200*xFactor,symY-0.068300,0,lc}; 
Point(11)={0.614300*xFactor,symY-0.050000,0,lc}; 
Point(12)={0.648500*xFactor,symY-0.068300,0,lc}; 
Point(13)={0.921500*xFactor,symY-0.068300,0,lc}; 
Point(14)={0.955600*xFactor,symY-0.034100,0,lc}; 
Point(15)={1.126300*xFactor,symY-0.034100,0,lc}; 
Point(16)={1.160400*xFactor,symY-0.170600,0,lc}; 
Point(17)={1.365200*xFactor,symY-0.170600,0,lc}; 
Point(18)={1.399300*xFactor,symY-0.204800,0,lc}; 
Point(19)={1.433400*xFactor,symY-0.238900,0,lc}; 
Point(20)={1.467600*xFactor,symY-0.273000,0,lc}; 
Point(21)={1.843000*xFactor,symY-0.273000,0,lc}; 
Point(22)={1.850200*xFactor,symY-0.265900,0,lc}; 
Point(23)={2.062200*xFactor,symY-0.409600,0,lc}; 
Point(24)={2.116000*xFactor,symY-0.409600,0,lc}; 
Point(25)={2.150200*xFactor,symY-0.443700,0,lc}; 
Point(26)={2.184300*xFactor,symY-0.443700,0,lc}; 
Point(27)={2.218400*xFactor,symY-0.511900,0,lc}; 
Point(28)={2.354900*xFactor,symY-0.511900,0,lc}; 
Point(29)={2.389100*xFactor,symY-0.546100,0,lc}; 
Point(30)={2.423200*xFactor,symY-0.546100,0,lc}; 
Point(31)={2.457300*xFactor,symY-0.580200,0,lc}; 
Point(32)={2.491500*xFactor,symY-0.580200,0,lc}; 
Point(33)={2.525600*xFactor,symY-0.614300,0,lc}; 
Point(34)={2.535500*xFactor,symY-0.614300,0,lc}; 
Point(35)={2.802500*xFactor,symY-0.750900,0,lc}; 
Point(36)={3.037500*xFactor,symY-0.750900,0,lc}; 
Point(37)={3.071700*xFactor,symY-0.716700,0,lc}; 
Point(38)={3.105800*xFactor,symY-0.716700,0,lc}; 
Point(39)={3.139900*xFactor,symY-0.750900,0,lc}; 
Point(40)={3.174100*xFactor,symY-0.750900,0,lc}; 
Point(41)={3.185100*xFactor,symY-0.761900,0,lc}; 
Point(42)={3.483300*xFactor,symY-0.785000,0,lc}; 
Point(43)={3.754300*xFactor,symY-0.785000,0,lc}; 
Point(44)={3.788400*xFactor,symY-0.750900,0,lc}; 
Point(45)={3.850200*xFactor,symY-0.750900,0,lc}; 
Point(46)={4.147500*xFactor,symY-0.716700,0,lc}; 
Point(47)={4.198000*xFactor,symY-0.716700,0,lc}; 
Point(48)={4.232100*xFactor,symY-0.682600,0,lc}; 
Point(49)={4.518800*xFactor,symY-0.682600,0,lc}; 
Point(50)={4.806800*xFactor,symY-0.711200,0,lc}; 
Point(51)={4.812300*xFactor,symY-0.716700,0,lc}; 
Point(52)={5.051200*xFactor,symY-0.716700,0,lc}; 
Point(53)={5.085300*xFactor,symY-0.682600,0,lc}; 
Point(54)={5.119500*xFactor,symY-0.648500,0,lc}; 
Point(55)={5.183300*xFactor,symY-0.648500,0,lc}; 
Point(56)={5.483300*xFactor,symY-0.648500,0,lc}; 
Point(57)={5.733800*xFactor,symY-0.648500,0,lc}; 
Point(58)={5.767900*xFactor,symY-0.682600,0,lc}; 
Point(59)={5.802000*xFactor,symY-0.785000,0,lc}; 
Point(60)={5.850000*xFactor,symY-0.785000,0,lc}; 
Point(61)={6.124500*xFactor,symY-0.701400,0,lc}; 
Point(62)={6.143300*xFactor,symY-0.682600,0,lc}; 
Point(63)={6.279900*xFactor,symY-0.682600,0,lc}; 
Point(64)={6.314000*xFactor,symY-0.648500,0,lc}; 
Point(65)={6.348100*xFactor,symY-0.648500,0,lc}; 
Point(66)={6.382300*xFactor,symY-0.614300,0,lc}; 
Point(67)={6.416400*xFactor,symY-0.626100,0,lc}; 
Point(68)={6.450500*xFactor,symY-0.682600,0,lc}; 
Point(69)={6.484600*xFactor,symY-0.716700,0,lc}; 
Point(70)={6.518800*xFactor,symY-0.750900,0,lc}; 
Point(71)={6.520300*xFactor,symY-0.752300,0,lc}; 
Point(72)={6.816700*xFactor,symY-0.785000,0,lc}; 
Point(73)={7.187300*xFactor,symY-0.785000,0,lc}; 
Point(74)={7.457200*xFactor,symY-0.666200,0,lc}; 
Point(75)={7.474400*xFactor,symY-0.580200,0,lc}; 
Point(76)={7.508500*xFactor,symY-0.580200,0,lc}; 
Point(77)={7.542700*xFactor,symY-0.546100,0,lc}; 
Point(78)={7.576800*xFactor,symY-0.511900,0,lc}; 
Point(79)={7.610900*xFactor,symY-0.511900,0,lc}; 
Point(80)={7.645100*xFactor,symY-0.409600,0,lc}; 
Point(81)={7.781600*xFactor,symY-0.409600,0,lc}; 
Point(82)={7.815700*xFactor,symY-0.375400,0,lc}; 
Point(83)={7.849800*xFactor,symY-0.307200,0,lc}; 
Point(84)={7.853200*xFactor,symY-0.303800,0,lc}; 
Point(85)={8.140700*xFactor,symY-0.221100,0,lc}; 
Point(86)={8.157000*xFactor,symY-0.204800,0,lc}; 
Point(87)={8.327600*xFactor,symY-0.204800,0,lc}; 
Point(88)={8.361800*xFactor,symY-0.170600,0,lc}; 
Point(89)={8.703100*xFactor,symY-0.170600,0,lc}; 
Point(90)={8.737200*xFactor,symY-0.136500,0,lc}; 
Point(91)={8.771300*xFactor,symY-0.136500,0,lc}; 
Point(92)={8.805500*xFactor,symY-0.102400,0,lc}; 
Point(93)={9.317400*xFactor,symY-0.102400,0,lc}; 
Point(94)={9.851500*xFactor,1.0,0,lc}; 
Point(95)={xMax*xFactor,1.0,0,lc}; 
Point(96)={2.000000*xFactor,symY-0.273000,0,lc}; 
Point(97)={2.666700*xFactor,symY-0.687100,0,lc}; 
Point(98)={3.333300*xFactor,symY-0.785000,0,lc}; 
Point(99)={4.000000*xFactor,symY-0.744000,0,lc}; 
Point(100)={4.666700*xFactor,symY-0.657600,0,lc}; 
Point(101)={5.333300*xFactor,symY-0.648500,0,lc}; 
Point(102)={6.000000*xFactor,symY-0.785000,0,lc}; 
Point(103)={6.666700*xFactor,symY-0.785000,0,lc}; 
Point(104)={7.333300*xFactor,symY-0.750900,0,lc}; 
Point(105)={8.000000*xFactor,symY-0.273000,0,lc}; 
Point(106)={1.915852*xFactor,symY-0.397212,0,lc}; 
Point(107)={2.598361*xFactor,symY-0.820676,0,lc}; 
Point(108)={3.321715*xFactor,symY-0.934547,0,lc}; 
Point(109)={4.017139*xFactor,symY-0.893000,0,lc}; 
Point(110)={4.651878*xFactor,symY-0.806867,0,lc}; 
Point(111)={5.333300*xFactor,symY-0.798500,0,lc}; 
Point(112)={6.043701*xFactor,symY-0.928474,0,lc}; 
Point(113)={6.650251*xFactor,symY-0.934099,0,lc}; 
Point(114)={7.393719*xFactor,symY-0.888201,0,lc}; 
Point(115)={8.041462*xFactor,symY-0.417137,0,lc}; 
Point(116)={0.000000*xFactor,0.00000,0,lc}; 
Point(117)={xMax*xFactor,0.00000,0,lc}; 
Line(1)={1,2};      
Line(2)={2,3};      
Line(3)={3,4};      
Line(4)={4,5};      
Line(5)={5,6};      
Line(6)={6,7};      
Line(7)={7,8};      
Line(8)={8,9};      
Line(9)={9,10};     
Line(10)={10,11};   
Line(11)={11,12};   
Line(12)={12,13};   
Line(13)={13,14};   
Line(14)={14,15};   
Line(15)={15,16};   
Line(16)={16,17};   
Line(17)={17,18};   
Line(18)={18,19};   
Line(19)={19,20};   
Line(20)={20,21};   
Line(21)={21,22};   
Line(22)={23,24};   
Line(23)={24,25};   
Line(24)={25,26};   
Line(25)={26,27};   
Line(26)={27,28};   
Line(27)={28,29};   
Line(28)={29,30}; 
Line(29)={30,31}; 
Line(30)={31,32}; 
Line(31)={32,33}; 
Line(32)={33,34}; 
Line(33)={35,36}; 
Line(34)={36,37}; 
Line(35)={37,38}; 
Line(36)={38,39}; 
Line(37)={39,40}; 
Line(38)={40,41}; 
Line(39)={42,43}; 
Line(40)={43,44}; 
Line(41)={44,45}; 
Line(42)={46,47}; 
Line(43)={47,48}; 
Line(44)={48,49}; 
Line(45)={50,51}; 
Line(46)={51,52}; 
Line(47)={52,53}; 
Line(48)={53,54}; 
Line(49)={54,55}; 
Line(50)={56,57}; 
Line(51)={57,58}; 
Line(52)={58,59}; 
Line(53)={59,60}; 
Line(54)={61,62}; 
Line(55)={62,63}; 
Line(56)={63,64}; 
Line(57)={64,65}; 
Line(58)={65,66}; 
Line(59)={66,67}; 
Line(60)={67,68}; 
Line(61)={68,69}; 
Line(62)={69,70}; 
Line(63)={70,71}; 
Line(64)={72,73}; 
Line(65)={74,75}; 
Line(66)={75,76}; 
Line(67)={76,77}; 
Line(68)={77,78}; 
Line(69)={78,79}; 
Line(70)={79,80}; 
Line(71)={80,81}; 
Line(72)={81,82}; 
Line(73)={82,83}; 
Line(74)={83,84}; 
Line(75)={85,86}; 
Line(76)={86,87}; 
Line(77)={87,88}; 
Line(78)={88,89}; 
Line(79)={89,90}; 
Line(80)={90,91}; 
Line(81)={91,92}; 
Line(82)={92,93}; 
Line(83)={93,94}; 
Line(84)={94,95}; 
Circle(85)={22,96,106}; 
Circle(86)={34,97,107}; 
Circle(87)={41,98,108}; 
Circle(88)={45,99,109}; 
Circle(89)={49,100,110}; 
Circle(90)={55,101,111}; 
Circle(91)={60,102,112}; 
Circle(92)={71,103,113}; 
Circle(93)={73,104,114}; 
Circle(94)={84,105,115}; 
Circle(95)={106,96,23}; 
Circle(96)={107,97,35}; 
Circle(97)={108,98,42}; 
Circle(98)={109,99,46}; 
Circle(99)={110,100,50}; 
Circle(100)={111,101,56}; 
Circle(101)={112,102,61}; 
Circle(102)={113,103,72}; 
Circle(103)={114,104,74}; 
Circle(104)={115,105,85}; 

Line(105)={116,117}; 
Line(106)={116,1}; 
Line(107)={117,95}; 
//+
Line Loop(108) = {105, 107, -84, -83, -82, -81, -80, -79, -78, -77, -76, -75, -104, -94, -74, -73, -72, -71, -70, -69, -68, -67, -66, -65, -103, -93, -64, -102, -92, -63, -62, -61, -60, -59, -58, -57, -56, -55, -54, -101, -91, -53, -52, -51, -50, -100, -90, -49, -48, -47, -46, -45, -99, -89, -44, -43, -42, -98, -88, -41, -40, -39, -97, -87, -38, -37, -36, -35, -34, -33, -96, -86, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -95, -85, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, -106};
//+
Plane Surface(109) = {108};

Physical Line("dirichlet1 noslip") = {85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84}; 
Physical Line("neumann1 outflow") = {107};
Physical Line("dirichlet1 inflow") = {106};
Physical Line("neumann1 symmetry") = {105};
Physical Line("dirichlet2 noslip") = {85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84}; 
Physical Line("neumann2 outflow") = {107};
Physical Line("dirichlet2 inflow") = {106};
Physical Line("neumann2 symmetry") = {105};
Physical Line("dirichlet3 bottom") = {105};
Physical Line("neumann3 outflow") = {107};
Physical Line("neumann3 inflow") = {106};
Physical Line("dirichlet3 top") = {85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84}; 
Physical Line("neumann4 bottom") = {105};
Physical Line("neumann4 outflow") = {107};
Physical Line("neumann4 inflow") = {106};
Physical Line("neumann4 top") = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84}; 
Physical Line("dirichlet4 stent") = {85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104}; 

Physical Surface(11) = {109};
