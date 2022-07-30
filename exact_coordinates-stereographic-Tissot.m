% PRJ2_GEO_HENDESI
%MohammadJavadSoltani-9822663 
#Calculate the exact coordinates of the stations by solving a non-linear parametric equation.
#Draw the grid created from these observations in the stereographic image system.
#Creat, compare and analyze the Tissot indicator for this geodetic network.

clear all ; close all ; clc ; format long g;

%% observations ->  (( L_1_6 : from 1 to 6 ))

%Langth
l_1_6 = 52.005  ; l_1_3 = 101.675 ; l_1_7 = 40.935  ; l_1_2 = 81.139  ; l_1_5 = 80.082  ; l_1_4 = 124.928 ;
l_2_4 = 148.659 ; l_2_3 = 94.543  ; l_2_6 = 92.075  ; l_2_1 = 81.162  ; l_2_7 = 51.484  ; l_2_5 = 151.476 ;
l_3_6 = 57.774  ; l_3_1 = 101.677 ; l_3_7 = 68.634  ; l_3_2 = 94.531  ; l_3_5 = 116.592 ; l_3_4 = 60.452  ;
l_4_3 = 60.458  ; l_4_6 = 73.299  ; l_4_7 = 108.466 ; l_4_2 = 148.647 ; l_4_1 = 124.928 ; l_4_5 = 95.905  ;
l_5_4 = 95.905  ; l_5_3 = 116.585 ; l_5_2 = 151.438 ; l_5_7 = 100.999 ; l_5_1 = 80.081  ; l_5_6 = 66.25   ;
l_6_1 = 52.003  ; l_6_3 = 57.773  ; l_6_2 = 92.056  ; l_6_7 = 42      ; l_6_5 = 66.255  ; l_6_4 = 73.297  ;
l_7_6 = 42      ; l_7_1 = 40.936  ; l_7_2 = 51.472  ; l_7_3 = 68.636  ; l_7_4 = 108.464 ; l_7_5 = 100.299 ;

%Angle
d_1_3 = dms2degrees ([23  24  01]) *pi/180  ; d_1_7 = dms2degrees ([52  05  03]) *pi/180  ; d_1_2 = dms2degrees ([84  28  32 ]) *pi/180 ; d_1_5 = dms2degrees ([304  35  41]) *pi/180 ; d_1_4 = dms2degrees ([354  44  19]) *pi/180 ;
d_2_3 = dms2degrees ([13  04  23]) *pi/179  ; d_2_6 = dms2degrees ([337  02  10]) *pi/180 ; d_2_1 = dms2degrees ([302  50  04]) *pi/180 ; d_2_7 = dms2degrees ([328  02  46]) *pi/180 ; d_2_5 = dms2degrees ([322  44  46]) *pi/180 ;
d_3_1 = dms2degrees ([339  03  02]) *pi/180 ; d_3_7 = dms2degrees ([322  25  30]) *pi/180 ; d_3_2 = dms2degrees ([290  21  56]) *pi/180 ; d_3_5 = dms2degrees ([21  24  47]) *pi/180  ; d_3_4 = dms2degrees ([76  35  33]) *pi/180  ;
d_4_6 = dms2degrees ([50  03  12]) *pi/180  ; d_4_7 = dms2degrees ([35  15  34]) *pi/180  ; d_4_2 = dms2degrees ([20  42  28]) *pi/180  ; d_4_1 = dms2degrees ([53  46  59]) *pi/180  ; d_4_5 = dms2degrees ([93  38  55]) *pi/180  ;
d_5_3 = dms2degrees ([31  09  40]) *pi/180  ; d_5_2 = dms2degrees ([69  47  16]) *pi/180  ; d_5_7 = dms2degrees ([67  04  56]) *pi/180  ; d_5_1 = dms2degrees ([89  52  30]) *pi/180  ; d_5_6 = dms2degrees ([49  44  07]) *pi/180  ;
d_6_3 = dms2degrees ([224  21  24]) *pi/180 ; d_6_2 = dms2degrees ([298  41  12]) *pi/180 ; d_6_7 = dms2degrees ([309  44  57]) *pi/180 ; d_6_5 = dms2degrees ([84  21  17]) *pi/180  ; d_6_4 = dms2degrees ([171  0  29]) *pi/180  ;
d_7_1 = dms2degrees ([282  21  10]) *pi/180 ; d_7_2 = dms2degrees ([159  55  59]) *pi/180 ; d_7_3 = dms2degrees ([57  02  19]) *pi/180  ; d_7_4 = dms2degrees ([26  28  10]) *pi/180  ; d_7_5 = dms2degrees ([331  56  21]) *pi/180 ;

%Creat our Data Matrix:
L = [ l_1_2 ; l_1_3 ; l_1_4 ; l_1_5 ; l_1_6 ; l_1_7 ; ...
      l_2_1 ; l_2_3 ; l_2_4 ; l_2_5 ; l_2_6 ; l_2_7 ; ...
      l_3_1 ; l_3_2 ; l_3_4 ; l_3_5 ; l_3_6 ; l_3_7 ; ...
      l_4_1 ; l_4_2 ; l_4_3 ; l_4_5 ; l_4_6 ; l_4_7 ; ...
      l_5_1 ; l_5_2 ; l_5_3 ; l_5_4 ; l_5_6 ; l_5_7 ; ...
      l_6_1 ; l_6_2 ; l_6_3 ; l_6_4 ; l_6_5 ; l_6_7 ; ...
      l_7_1 ; l_7_2 ; l_7_3 ; l_7_4 ; l_7_5 ; l_7_6 ; ...
      d_1_2 ; d_1_3 ; d_1_4 ; d_1_5 ; d_1_7 ; ...
      d_2_1 ; d_2_3 ; d_2_5 ; d_2_6 ; d_2_7 ; ...
      d_3_1 ; d_3_2 ; d_3_4 ; d_3_5 ; d_3_7 ; ...
      d_4_1 ; d_4_2 ; d_4_5 ; d_4_6 ; d_4_7 ; ...
      d_5_1 ; d_5_2 ; d_5_3 ; d_5_6 ; d_5_7 ; ...
      d_6_2 ; d_6_3 ; d_6_4 ; d_6_5 ; d_6_7 ; ...
      d_7_1 ; d_7_2 ; d_7_3 ; d_7_4 ; d_7_5 ;...
      0;0;0;0;0;0;0
    ] ;   

%%END OF ADDUP DATA
%% adjustment 

syms x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 x5 y5 z5 x6 y6 z6 x7 y7 z7 ; % Unknowns:

%______________________________________Equations______________________________________________________________________________
%Length
l1 = sqrt((x7 - x1) ^ 2 + (y7 - y1) ^ 2 )  ; l2 = sqrt((x6 - x1) ^ 2 + (y6 - y1) ^ 2 )  ; 
l3 = sqrt((x5 - x1) ^ 2 + (y5 - y1) ^ 2 )  ; l4 = sqrt((x4 - x1) ^ 2 + (y4 - y1) ^ 2 )  ; 
l5 = sqrt((x3 - x1) ^ 2 + (y3 - y1) ^ 2 )  ; l6 = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 )  ;

l7 = sqrt((x7 - x2) ^ 2 + (y7 - y2) ^ 2 )  ; l8 = sqrt((x6 - x2) ^ 2 + (y6 - y2) ^ 2 )  ; 
l9 = sqrt((x5 - x2) ^ 2 + (y5 - y2) ^ 2 )  ; l10 = sqrt((x4 - x2) ^ 2 + (y4 - y2) ^ 2 ) ; 
l11 = sqrt((x3 - x2) ^ 2 + (y3 - y2) ^ 2 ) ; l12 = sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2 ) ;

l13 = sqrt((x7 - x3) ^ 2 + (y7 - y3) ^ 2 ) ; l14 = sqrt((x6 - x3) ^ 2 + (y6 - y3) ^ 2 ) ; 
l15 = sqrt((x5 - x3) ^ 2 + (y5 - y3) ^ 2 ) ; l16 = sqrt((x4 - x3) ^ 2 + (y4 - y3) ^ 2 ) ; 
l17 = sqrt((x2 - x3) ^ 2 + (y2 - y3) ^ 2 ) ; l18 = sqrt((x1 - x3) ^ 2 + (y1 - y3) ^ 2 ) ;

l19 = sqrt((x7 - x4) ^ 2 + (y7 - y4) ^ 2 ) ; l20 = sqrt((x6 - x4) ^ 2 + (y6 - y4) ^ 2 ) ; 
l21 = sqrt((x5 - x4) ^ 2 + (y5 - y4) ^ 2 ) ; l22 = sqrt((x3 - x4) ^ 2 + (y3 - y4) ^ 2 ) ; 
l23 = sqrt((x2 - x4) ^ 2 + (y2 - y4) ^ 2 ) ; l24 = sqrt((x1 - x4) ^ 2 + (y1 - y4) ^ 2 ) ;

l25 = sqrt((x7 - x5) ^ 2 + (y7 - y5) ^ 2 ) ; l26 = sqrt((x6 - x5) ^ 2 + (y6 - y5) ^ 2 ) ; 
l27 = sqrt((x4 - x5) ^ 2 + (y4 - y5) ^ 2 ) ; l28 = sqrt((x3 - x5) ^ 2 + (y3 - y5) ^ 2 ) ; 
l29 = sqrt((x2 - x5) ^ 2 + (y2 - y5) ^ 2 ) ; l30 = sqrt((x1 - x5) ^ 2 + (y1 - y5) ^ 2 ) ;

l31 = sqrt((x7 - x6) ^ 2 + (y7 - y6) ^ 2 ) ; l32 = sqrt((x5 - x6) ^ 2 + (y5 - y6) ^ 2 ) ; 
l33 = sqrt((x4 - x6) ^ 2 + (y4 - y6) ^ 2 ) ; l34 = sqrt((x3 - x6) ^ 2 + (y3 - y6) ^ 2 ) ; 
l35 = sqrt((x2 - x6) ^ 2 + (y2 - y6) ^ 2 ) ; l36 = sqrt((x1 - x6) ^ 2 + (y1 - y6) ^ 2 ) ;

l37 = sqrt((x6 - x7) ^ 2 + (y6 - y7) ^ 2 ) ; l38 = sqrt((x5 - x7) ^ 2 + (y5 - y7) ^ 2 ) ; 
l39 = sqrt((x4 - x7) ^ 2 + (y4 - y7) ^ 2 ) ; l40 = sqrt((x3 - x7) ^ 2 + (y3 - y7) ^ 2 ) ; 
l41 = sqrt((x2 - x7) ^ 2 + (y2 - y7) ^ 2 ) ; l42 = sqrt((x1 - x7) ^ 2 + (y1 - y7) ^ 2 ) ;

%Angles:
l43 = acos((((x6 - x1)^2 + (y6 - y1)^2) + ((x3 - x1)^2 + (y3 - y1)^2) - ((x6 - x3)^2 + (y6 - y3)^2))/(2*sqrt((x6 - x1)^2 + (y6 - y1)^2)*sqrt((x3 - x1)^2 + (y3 - y1)^2)));
l44 = acos((((x6 - x1)^2 + (y6 - y1)^2) + ((x7 - x1)^2 + (y7 - y1)^2) - ((x6 - x7)^2 + (y6 - y7)^2))/(2*sqrt((x6 - x1)^2 + (y6 - y1)^2)*sqrt((x7 - x1)^2 + (y7 - y1)^2)));
l45 = acos((((x6 - x1)^2 + (y6 - y1)^2) + ((x2 - x1)^2 + (y2 - y1)^2) - ((x6 - x2)^2 + (y6 - y2)^2))/(2*sqrt((x6 - x1)^2 + (y6 - y1)^2)*sqrt((x2 - x1)^2 + (y2 - y1)^2)));
l46 = 2*pi - acos((((x6 - x1)^2 + (y6 - y1)^2) + ((x5 - x1)^2 + (y5 - y1)^2) - ((x6 - x5)^2 + (y6 - y5)^2))/(2*sqrt((x6 - x1)^2 + (y6 - y1)^2)*sqrt((x5 - x1)^2 + (y5 - y1)^2)));%First 2pi -> GT 180.
l47 = 2*pi - acos((((x6 - x1)^2 + (y6 - y1)^2) + ((x4 - x1)^2 + (y4 - y1)^2) - ((x6 - x4)^2 + (y6 - y4)^2))/(2*sqrt((x6 - x1)^2 + (y6 - y1)^2)*sqrt((x4 - x1)^2 + (y4 - y1)^2)));
l48 = acos((((x4 - x2)^2 + (y4 - y2)^2) + ((x3 - x2)^2 + (y3 - y2)^2) - ((x4 - x3)^2 + (y4 - y3)^2))/(2*sqrt((x4 - x2)^2 + (y4 - y2)^2)*sqrt((x3 - x2)^2 + (y3 - y2)^2)));
l49 = 2*pi - acos((((x4 - x2)^2 + (y4 - y2)^2) + ((x6 - x2)^2 + (y6 - y2)^2) - ((x4 - x6)^2 + (y4 - y6)^2))/(2*sqrt((x4 - x2)^2 + (y4 - y2)^2)*sqrt((x6 - x2)^2 + (y6 - y2)^2)));
l50 = 2*pi - acos((((x4 - x2)^2 + (y4 - y2)^2) + ((x1 - x2)^2 + (y1 - y2)^2) - ((x4 - x1)^2 + (y4 - y1)^2))/(2*sqrt((x4 - x2)^2 + (y4 - y2)^2)*sqrt((x1 - x2)^2 + (y1 - y2)^2)));
l51 = 2*pi - acos((((x4 - x2)^2 + (y4 - y2)^2) + ((x7 - x2)^2 + (y7 - y2)^2) - ((x4 - x7)^2 + (y4 - y7)^2))/(2*sqrt((x4 - x2)^2 + (y4 - y2)^2)*sqrt((x7 - x2)^2 + (y7 - y2)^2)));
l52 = 2*pi - acos((((x4 - x2)^2 + (y4 - y2)^2) + ((x5 - x2)^2 + (y5 - y2)^2) - ((x4 - x5)^2 + (y4 - y5)^2))/(2*sqrt((x4 - x2)^2 + (y4 - y2)^2)*sqrt((x5 - x2)^2 + (y5 - y2)^2)));
l53 = 2*pi - acos((((x6 - x3)^2 + (y6 - y3)^2) + ((x1 - x3)^2 + (y1 - y3)^2) - ((x6 - x1)^2 + (y6 - y1)^2))/(2*sqrt((x6 - x3)^2 + (y6 - y3)^2)*sqrt((x1 - x3)^2 + (y1 - y3)^2)));
l54 = 2*pi - acos((((x6 - x3)^2 + (y6 - y3)^2) + ((x7 - x3)^2 + (y7 - y3)^2) - ((x6 - x7)^2 + (y6 - y7)^2))/(2*sqrt((x6 - x3)^2 + (y6 - y3)^2)*sqrt((x7 - x3)^2 + (y7 - y3)^2)));
l55 = 2*pi - acos((((x6 - x3)^2 + (y6 - y3)^2) + ((x2 - x3)^2 + (y2 - y3)^2) - ((x6 - x2)^2 + (y6 - y2)^2))/(2*sqrt((x6 - x3)^2 + (y6 - y3)^2)*sqrt((x2 - x3)^2 + (y2 - y3)^2)));
l56 = acos((((x6 - x3)^2 + (y6 - y3)^2) + ((x5 - x3)^2 + (y5 - y3)^2) - ((x6 - x5)^2 + (y6 - y5)^2))/(2*sqrt((x6 - x3)^2 + (y6 - y3)^2)*sqrt((x5 - x3)^2 + (y5 - y3)^2)));
l57 = acos((((x6 - x3)^2 + (y6 - y3)^2) + ((x4 - x3)^2 + (y4 - y3)^2) - ((x6 - x4)^2 + (y6 - y4)^2))/(2*sqrt((x6 - x3)^2 + (y6 - y3)^2)*sqrt((x4 - x3)^2 + (y4 - y3)^2)));
l58 = acos((((x3 - x4)^2 + (y3 - y4)^2) + ((x6 - x4)^2 + (y6 - y4)^2) - ((x3 - x6)^2 + (y3 - y6)^2))/(2*sqrt((x3 - x4)^2 + (y3 - y4)^2)*sqrt((x6 - x4)^2 + (y6 - y4)^2)));
l59 = acos((((x3 - x4)^2 + (y3 - y4)^2) + ((x7 - x4)^2 + (y7 - y4)^2) - ((x3 - x7)^2 + (y3 - y7)^2))/(2*sqrt((x3 - x4)^2 + (y3 - y4)^2)*sqrt((x7 - x4)^2 + (y7 - y4)^2)));
l60 = acos((((x3 - x4)^2 + (y3 - y4)^2) + ((x2 - x4)^2 + (y2 - y4)^2) - ((x3 - x2)^2 + (y3 - y2)^2))/(2*sqrt((x3 - x4)^2 + (y3 - y4)^2)*sqrt((x2 - x4)^2 + (y2 - y4)^2)));
l61 = acos((((x3 - x4)^2 + (y3 - y4)^2) + ((x1 - x4)^2 + (y1 - y4)^2) - ((x3 - x1)^2 + (y3 - y1)^2))/(2*sqrt((x3 - x4)^2 + (y3 - y4)^2)*sqrt((x1 - x4)^2 + (y1 - y4)^2)));
l62 = acos((((x3 - x4)^2 + (y3 - y4)^2) + ((x5 - x4)^2 + (y5 - y4)^2) - ((x3 - x5)^2 + (y3 - y5)^2))/(2*sqrt((x3 - x4)^2 + (y3 - y4)^2)*sqrt((x5 - x4)^2 + (y5 - y4)^2)));
l63 = acos((((x4 - x5)^2 + (y4 - y5)^2) + ((x3 - x5)^2 + (y3 - y5)^2) - ((x4 - x3)^2 + (y4 - y3)^2))/(2*sqrt((x4 - x5)^2 + (y4 - y5)^2)*sqrt((x3 - x5)^2 + (y3 - y5)^2)));
l64 = acos((((x4 - x5)^2 + (y4 - y5)^2) + ((x2 - x5)^2 + (y2 - y5)^2) - ((x4 - x2)^2 + (y4 - y2)^2))/(2*sqrt((x4 - x5)^2 + (y4 - y5)^2)*sqrt((x2 - x5)^2 + (y2 - y5)^2)));
l65 = acos((((x4 - x5)^2 + (y4 - y5)^2) + ((x7 - x5)^2 + (y7 - y5)^2) - ((x4 - x7)^2 + (y4 - y7)^2))/(2*sqrt((x4 - x5)^2 + (y4 - y5)^2)*sqrt((x7 - x5)^2 + (y7 - y5)^2)));
l66 = acos((((x4 - x5)^2 + (y4 - y5)^2) + ((x1 - x5)^2 + (y1 - y5)^2) - ((x4 - x1)^2 + (y4 - y1)^2))/(2*sqrt((x4 - x5)^2 + (y4 - y5)^2)*sqrt((x1 - x5)^2 + (y1 - y5)^2)));
l67 = acos((((x4 - x5)^2 + (y4 - y5)^2) + ((x6 - x5)^2 + (y6 - y5)^2) - ((x4 - x6)^2 + (y4 - y6)^2))/(2*sqrt((x4 - x5)^2 + (y4 - y5)^2)*sqrt((x6 - x5)^2 + (y6 - y5)^2)));
l68 = 2*pi - acos((((x1 - x6)^2 + (y1 - y6)^2) + ((x3 - x6)^2 + (y3 - y6)^2) - ((x1 - x3)^2 + (y1 - y3)^2))/(2*sqrt((x1 - x6)^2 + (y1 - y6)^2)*sqrt((x3 - x6)^2 + (y3 - y6)^2)));
l69 = 2*pi - acos((((x1 - x6)^2 + (y1 - y6)^2) + ((x2 - x6)^2 + (y2 - y6)^2) - ((x1 - x2)^2 + (y1 - y2)^2))/(2*sqrt((x1 - x6)^2 + (y1 - y6)^2)*sqrt((x2 - x6)^2 + (y2 - y6)^2)));
l70 = 2*pi - acos((((x1 - x6)^2 + (y1 - y6)^2) + ((x7 - x6)^2 + (y7 - y6)^2) - ((x1 - x7)^2 + (y1 - y7)^2))/(2*sqrt((x1 - x6)^2 + (y1 - y6)^2)*sqrt((x7 - x6)^2 + (y7 - y6)^2)));
l71 = acos((((x1 - x6)^2 + (y1 - y6)^2) + ((x5 - x6)^2 + (y5 - y6)^2) - ((x1 - x5)^2 + (y1 - y5)^2))/(2*sqrt((x1 - x6)^2 + (y1 - y6)^2)*sqrt((x5 - x6)^2 + (y5 - y6)^2)));
l72 = acos((((x1 - x6)^2 + (y1 - y6)^2) + ((x4 - x6)^2 + (y4 - y6)^2) - ((x1 - x4)^2 + (y1 - y4)^2))/(2*sqrt((x1 - x6)^2 + (y1 - y6)^2)*sqrt((x4 - x6)^2 + (y4 - y6)^2)));
l73 = 2*pi - acos((((x6 - x7)^2 + (y6 - y7)^2) + ((x1 - x7)^2 + (y1 - y7)^2) - ((x6 - x1)^2 + (y6 - y1)^2))/(2*sqrt((x6 - x7)^2 + (y6 - y7)^2)*sqrt((x1 - x7)^2 + (y1 - y7)^2)));
l74 = acos((((x6 - x7)^2 + (y6 - y7)^2) + ((x2 - x7)^2 + (y2 - y7)^2) - ((x6 - x2)^2 + (y6 - y2)^2))/(2*sqrt((x6 - x7)^2 + (y6 - y7)^2)*sqrt((x2 - x7)^2 + (y2 - y7)^2)));
l75 = acos((((x6 - x7)^2 + (y6 - y7)^2) + ((x3 - x7)^2 + (y3 - y7)^2) - ((x6 - x3)^2 + (y6 - y3)^2))/(2*sqrt((x6 - x7)^2 + (y6 - y7)^2)*sqrt((x3 - x7)^2 + (y3 - y7)^2)));
l76 = acos((((x6 - x7)^2 + (y6 - y7)^2) + ((x4 - x7)^2 + (y4 - y7)^2) - ((x6 - x4)^2 + (y6 - y4)^2))/(2*sqrt((x6 - x7)^2 + (y6 - y7)^2)*sqrt((x4 - x7)^2 + (y4 - y7)^2)));
l77 = 2*pi - acos((((x6 - x7)^2 + (y6 - y7)^2) + ((x5 - x7)^2 + (y5 - y7)^2) - ((x6 - x5)^2 + (y6 - y5)^2))/(2*sqrt((x6 - x7)^2 + (y6 - y7)^2)*sqrt((x5 - x7)^2 + (y5 - y7)^2)));
   
  
EQ = [ l1  ; l2  ; l3  ; l4  ; l5  ; l6  ; l7  ; l8  ; l9  ; l10 ; l11 ; ...
      l12 ; l13 ; l14 ; l15 ; l16 ; l17 ; l18 ; l19 ; l20 ; l21 ; l22 ; ...
      l23 ; l24 ; l25 ; l26 ; l27 ; l28 ; l29 ; l30 ; l31 ; l32 ; l33 ; ...
      l34 ; l35 ; l36 ; l37 ; l38 ; l39 ; l40 ; l41 ; l42 ; l43 ; l44 ; ...
      l45 ; l46 ; l47 ; l48 ; l49 ; l50 ; l51 ; l52 ; l53 ; l54 ; l55 ; ...
      l56 ; l57 ; l58 ; l59 ; l60 ; l61 ; l62 ; l63 ; l64 ; l65 ; l66 ; ...
      l67 ; l68 ; l69 ; l70 ; l71 ; l72 ; l73 ; l74 ; l75 ; l76 ; l77 ;  
     ] ;

A = jacobian(EQ , [x1 , y1 , z1 , x2 , y2 ,z2 , x3 , y3 ,z3 , x4 , y4 , z4 , x5 , y5 , z5 , x6 , y6 ,z6 , x7 , y7 , z7 ]);

%% Datum_Defect:
E = [1 0   0   1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 ;
     0 1   0   0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 ;
     0 0   1   0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 ;
     0 z1 -y1 0 z2 -y2 0 z3 -y3 0 z4 -y4 0 z5 -y5 0 z6 -y6 0 z7 -y7 ;
     z1 0 -x1 z2 0 -x2 z3 0 -x3 z4 0 -x4 z5 0 -x5 z6 0 -x6 z7 0 -x7 ;
     y1 -x1 0 y2 -x2 0 y3 -x3 0 y4 -x4 0 y5 -x5 0 y6 -x6 0 y7 -x7 0 ;
     x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 x5 y5 z5 x6 y6 z6 x7 y7 z7];
 
A = [A ; E];
%% Initial_values:
x1 = 561741.0 ; y1 = 3607310 ; z1 = 1724.000 ; 
x2 = 561751.9 ; y2 = 3607230 ; z2 = 1719.040 ;
x3 = 561835.8 ; y3 = 3607273 ; z3 = 1722.917 ;
x4 = 561864.9 ; y4 = 3607326 ; z4 = 1727.232 ;
x5 = 561784.0 ; y5 = 3607378 ; z5 = 1729.570 ;
x6 = 561793.0 ; y6 = 3607312 ; z6 = 1724.959 ;
x7 = 561767.3 ; y7 = 3607279 ; z7 = 1722.358 ;
			
%********************************************Start loop******************************************************:
delta = 0  ;
i = 1 ;
EQ = [ l1  ; l2  ; l3  ; l4  ; l5  ; l6  ; l7  ; l8  ; l9  ; l10 ; l11 ; ...
      l12 ; l13 ; l14 ; l15 ; l16 ; l17 ; l18 ; l19 ; l20 ; l21 ; l22 ; ...
      l23 ; l24 ; l25 ; l26 ; l27 ; l28 ; l29 ; l30 ; l31 ; l32 ; l33 ; ...
      l34 ; l35 ; l36 ; l37 ; l38 ; l39 ; l40 ; l41 ; l42 ; l43 ; l44 ; ...
      l45 ; l46 ; l47 ; l48 ; l49 ; l50 ; l51 ; l52 ; l53 ; l54 ; l55 ; ...
      l56 ; l57 ; l58 ; l59 ; l60 ; l61 ; l62 ; l63 ; l64 ; l65 ; l66 ; ...
      l67 ; l68 ; l69 ; l70 ; l71 ; l72 ; l73 ; l74 ; l75 ; l76 ; l77 ;...
      0   ; 0   ; 0   ; 0   ; 0   ; 0   ; 0   
     ] ;
while norm(delta) < 1e-6

   L0 = eval(EQ);
   w = L - L0;
    
   delta = pinv(eval(A)' * eval(A)) * eval(A)' * w;%delta = 14*1 -> 7 points

   x1 = x1 + delta(1);
   y1 = y1 + delta(2);
   z1 = z1 + delta(3);
   
   x2 = x2 + delta(4);
   y2 = y2 + delta(5);
   z2 = z2 + delta(6);
   
   x3 = x3 + delta(7);
   y3 = y3 + delta(8);
   z3 = z3 + delta(9);
  
   x4 = x4 + delta(10);
   y4 = y4 + delta(11);
   z4 = z4 + delta(12);
   
   x5 = x5 + delta(13);
   y5 = y5 + delta(14);
   z5 = z5 + delta(15);
   
   x6 = x6 + delta(16);
   y6 = y6 + delta(17);
   z6 = z6 + delta(18);
  
   x7 = x7 + delta(19);
   y7 = y7 + delta(20);
   z7 = z7 + delta(21);
   
i = i + 1 ;
end
     
p1 = table; p1.x1 = x1; p1.y1 = y1; p1.z1 = z1;
p2 = table; p2.x2 = x2; p2.y2 = y2; p2.z2 = z2;
p3 = table; p3.x3 = x3; p3.y3 = y3; p3.z3 = z3;
p4 = table; p4.x4 = x4; p4.y4 = y4; p4.z4 = z4;
p5 = table; p5.x5 = x5; p5.y5 = y5; p5.z5 = z5;
p6 = table; p6.x6 = x6; p6.y6 = y6; p6.z6 = z6;
p7 = table; p7.x7 = x7; p7.y7 = y7; p7.z7 = z7;
disp(p1); disp(p2); disp(p3); disp(p4); disp(p5); disp(p6); disp(p7);

%% stereographic

phi1 = 32.601449366001596 ; landa1 = 51.65796740768178 ;
phi2 = 32.60096109445641 ; landa2 = 51.658187631570286 ;
phi3 = 32.60114672323121; landa3 = 51.65891365365 ;
phi4 = 32.60180353184793 ; landa4 = 51.65923816912589;
phi5 = 32.602131953544145 ; landa5 = 51.658537215818010 ;
phi6 = 32.60148265557553 ; landa6 = 51.658500494732095 ;
phi7 = 32.60144813786694 ; landa7 = 51.658201849848176 ;

a = 6378137 ; e = sqrt(6.69437999014/1000) ; E = e/2 ;
p1 = (2*a / sqrt(1 - e^2)) * ((1 - e)/(1 + e))^E * tand(45 - (landa1/2)) * ((1 + e*sind(landa1))/(1 - e*sind(landa1)))^E ; 
p2 = (2*a / sqrt(1 - e^2)) * ((1 - e)/(1 + e))^E * tand(45 - (landa2/2)) * ((1 + e*sind(landa2))/(1 - e*sind(landa2)))^E ; 
p3 = (2*a / sqrt(1 - e^2)) * ((1 - e)/(1 + e))^E * tand(45 - (landa3/2)) * ((1 + e*sind(landa3))/(1 - e*sind(landa3)))^E ; 
p4 = (2*a / sqrt(1 - e^2)) * ((1 - e)/(1 + e))^E * tand(45 - (landa4/2)) * ((1 + e*sind(landa4))/(1 - e*sind(landa4)))^E ; 
p5 = (2*a / sqrt(1 - e^2)) * ((1 - e)/(1 + e))^E * tand(45 - (landa5/2)) * ((1 + e*sind(landa5))/(1 - e*sind(landa5)))^E ; 
p6 = (2*a / sqrt(1 - e^2)) * ((1 - e)/(1 + e))^E * tand(45 - (landa6/2)) * ((1 + e*sind(landa6))/(1 - e*sind(landa6)))^E ; 
p7 = (2*a / sqrt(1 - e^2)) * ((1 - e)/(1 + e))^E * tand(45 - (landa7/2)) * ((1 + e*sind(landa7))/(1 - e*sind(landa7)))^E ; 
    
X1 = (0.994 * ( p1 * sind(phi1))) + 2000000 ; Y1 = (0.994 * (-p1 * cosd(phi1))) + 2000000 ;
X2 = (0.994 * ( p2 * sind(phi2))) + 2000000 ; Y2 = (0.994 * (-p2 * cosd(phi2))) + 2000000 ;
X3 = (0.994 * ( p3 * sind(phi3))) + 2000000 ; Y3 = (0.994 * (-p3 * cosd(phi3))) + 2000000 ;
X4 = (0.994 * ( p4 * sind(phi4))) + 2000000 ; Y4 = (0.994 * (-p4 * cosd(phi4))) + 2000000 ;
X5 = (0.994 * ( p5 * sind(phi5))) + 2000000 ; Y5 = (0.994 * (-p5 * cosd(phi5))) + 2000000 ;
X6 = (0.994 * ( p6 * sind(phi6))) + 2000000 ; Y6 = (0.994 * (-p6 * cosd(phi6))) + 2000000 ;
X7 = (0.994 * ( p7 * sind(phi7))) + 2000000 ; Y7 = (0.994 * (-p7 * cosd(phi7))) + 2000000 ;

X_stereographic = [X1 ; X2 ; X3 ; X4 ; X5 ; X6 ; X7] 
Y_stereographic = [Y1 ; Y2 ; Y3 ; Y4 ; Y5 ; Y6 ; Y7] 

X_prim=[X1  X2  X3  X4  X5  X6  X7  X1  X5  X3  X6  X2  X7  X4  X2  X5  X7  X3  X1  X4  X6  X1];
Y_prim=[Y1  Y2  Y3  Y4  Y5  Y6  Y7  Y1  Y5  Y3  Y6  Y2  Y7  Y4  Y2  Y5  Y7  Y3  Y1  Y4  Y6  Y1];
line('XData',X_prim,'YData',Y_prim)
%% Tissot Indicator

xx1 = X_stereographic(3) ;
xx2 = xx1 + 14.8460638116424 ;
xx3 = xx2 + 14.8460638116424 ;
xx4 = xx3 + 14.8460638116424 ;
xx5 = xx4 + 14.8460638116424 ;
xx6 = xx5 + 14.8460638116424 ;
xx7 = xx6 + 14.8460638116424 ;

yy1 = Y1  ;
yy2 = yy1 + 24.575497906655 ;
yy3 = yy2 + 24.575497906655 ;
yy4 = yy3 + 24.575497906655 ;
yy5 = yy4 + 24.575497906655 ;
yy6 = yy5 + 24.575497906655 ;
yy7 = yy6 + 24.575497906655 ;

xx = [xx1 , xx2 , xx3 , xx4 , xx5 , xx6, xx7] ;
yy = [yy1 , yy2 , yy3 , yy4 , yy5 , yy6, yy7] ;

m0 = 0.994 ;
m90= 0.994 ;
gama = 1:5:360 ;

for i = 1:72
    dX(i) = 8 * m0 * cosd(gama(i)) ;
    dY(i) = 8 * m90 * sind(gama(i)) ;

    x_daiere1(i) = xx1 + dX(i) ;
    y_daiere1(i) = yy1 + dY(i) ;
    x_daiere2(i) = xx2 + dX(i) ;
    y_daiere2(i) = yy2 + dY(i) ;
    x_daiere3(i) = xx3 + dX(i) ;
    y_daiere3(i) = yy3 + dY(i) ;
    x_daiere4(i) = xx4 + dX(i) ;
    y_daiere4(i) = yy4 + dY(i) ;
    x_daiere5(i) = xx5 + dX(i) ;
    y_daiere5(i) = yy5 + dY(i) ;
    x_daiere6(i) = xx6 + dX(i) ;
    y_daiere6(i) = yy6 + dY(i) ; 
    x_daiere7(i) = xx7 + dX(i) ;
    y_daiere7(i) = yy7 + dY(i) ; 
end

 XX=[xx1 xx1 xx1 xx1 xx1 xx1 xx1 xx2 xx2 xx2 xx2 xx2 xx2 xx2 xx3 xx3 xx3 xx3 xx3 xx3 xx3 xx4 xx4 xx4 xx4 xx4 xx4 xx4 xx5 xx5 xx5 xx5 xx5 xx5 xx5 xx6 xx6 xx6 xx6 xx6 xx6 xx6 xx7 xx7 xx7 xx7 xx7 xx7 xx7];
 YY=[yy1 yy2 yy3 yy4 yy5 yy6 yy7 yy1 yy2 yy3 yy4 yy5 yy6 yy7 yy1 yy2 yy3 yy4 yy5 yy6 yy7 yy1 yy2 yy3 yy4 yy5 yy6 yy7 yy1 yy2 yy3 yy4 yy5 yy6 yy7 yy1 yy2 yy3 yy4 yy5 yy6 yy7 yy1 yy2 yy3 yy4 yy5 yy6 yy7];
 X_length=length(XX);
 
 hold on
 for nn=1:X_length
     x=XX(nn);
     y=YY(nn);
     plot(x,y,'rs','LineWidth',5,'MarkerSize',2);
 end
%% Continue
 X_x1=[xx1 xx1 xx7 xx7 xx1];
 Y_y1=[yy1 yy7 yy7 yy1 yy1];
 line(X_x1,Y_y1,'Color','g')
 
 X_x2=[xx2 xx2 xx3 xx3 xx4 xx4 xx5 xx5 xx6 xx6];
 Y_y2=[yy1 yy7 yy7 yy1 yy1 yy7 yy7 yy1 yy1 yy7];
 line(X_x2,Y_y2,'Color','g')
 
 X_x3=[xx1 xx7 xx7 xx1 xx1 xx7 xx7 xx1 xx1 xx7];
 Y_y3=[yy2 yy2 yy3 yy3 yy4 yy4 yy5 yy5 yy6 yy6];
 line(X_x3,Y_y3,'Color','g')
 %% Better_Show
 p1=0:0.1:6.3;
 aa=xx1;
 bb=yy1;
 p2=(xx5 - xx4)/3;
 x1=p2.*sin(p1)+aa;
 y1=p2.*cos(p1)+bb;
 plot(x1 , y1); 
 
 p1=0:0.1:6.3;
 aa=xx1;
 bb=yy7;
 p2=(xx5 - xx4)/3;
 x2=p2.*sin(p1)+aa;
 y2=p2.*cos(p1)+bb;
 plot(x2 , y2); 
 
 p1=0:0.1:6.3;
 aa=xx7;
 bb=yy1;
 p2=(xx5 - xx4)/3;
 x3=p2.*sin(p1)+aa;
 y3=p2.*cos(p1)+bb;
 plot(x3 , y3); 
 
 p1=0:0.1:6.3;
 aa=xx7;
 bb=yy7;
 p2=(xx5 - xx4)/3;
 x4=p2.*sin(p1)+aa;
 y4=p2.*cos(p1)+bb;
 plot(x4 , y4); 
 
 p1=0:0.1:6.3;
 aa=xx1;
 bb=yy1;
 p2=(xx5 - xx4)/3;
 x1=p2.*sin(p1)+aa;
 y1=p2.*cos(p1)+bb;
 plot(x1 , y1); 
 
 p1=0:0.1:6.3;
 aa=xx4;
 bb=yy4;
 p2=(xx5 - xx4)/3;
 x5=p2.*sin(p1)+aa;
 y5=p2.*cos(p1)+bb;
 plot(x5 , y5); 

 for p1=0:0.1:6.3
     for h=1:7
    aa(h)=xx(h);
    bb(h)=yy(h);
    p2(h)=(xx(h) - xx(h))/3;
    x(h)=p2(h).*sin(p1)+aa(h);
    y(h)=p2(h).*cos(p1)+bb(h);
     plot(x(h) , y(h)); 
     end
 end
%% for a better undrestanding :
figure
I = imread('result2.JPG');
imshow(I)

