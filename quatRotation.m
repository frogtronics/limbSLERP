(* ::Package:: *)

BeginPackage["quatRotation`",{"Quaternions`"}]

qHelp::usage ="type this to get help and example code"

qCon::usage = "[q]  Find the conjugate of a quaternion"

qMult::usage = "[q1, q2]  Multiply two quaternions together"

qRot::usage = "[q, 3vector]  Rotates a given 3D vector using a quaternion"

qqVec::usage= "[q1, q2] Computes the quaternion to rotate one 3d vector onto another"

qAngle::usage = "[q] Extracts the scalar angle from a unit quaternion. The 3D orientation of
an object can be represented as an axis of rotation and an angle about that axis"

qAxis::usage = "[q] Extracts the unit vector rotation axis of a unit quaternion.  The 3D orientation of
an object can be represented as an axis of rotation and an angle about that axis"

qtoR::usage = "[q] If you really have to do it, this will give you a rotation matrix from a quaternion.
But when you convert, you lose the benefits of quaternions"

qtoEuler::usage = "[q, 3vector]. Output a series of Euler angles from a quaternion which
is painful. q is the quaternion and 3vector is the rotation indices determining rotation order.  For example,
x-y-z would be {1, 2, 3}, x-z-x would be {1,3,1}, etc. This function uses rather unconventional 
math (ref. 1).  For more conventional formulae see ref. 2. Note: 321 sequence doesn't appear to work, so 
use qtoEuler321 instead. Have fun."

qtoEuler321::usage = "[q]. 3-2-1 rotation sequence from quat. 
Based on wikipedia which seems to work for this sequence"

qDot::usage = "[axisUnitVector3, order=3] Time derivative of a quaternion.  Ultimately useful for
calculating 3D rotational velocities without having to use messy Euler angles.  axisUnitVector3 is a 3D
unit vector for the quaternion's rotation axis.  Note: theta must be kept as a symbol!"


qHelp::usage = "Print example code"

(*Print["Code by Chris Richards, Royal Veterinary College, 2015 

Dependencies:  Quaternions package
References: 

ref 1.  http://noelhughes.net/Quaternion_to_DCM.html

ref 2. https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions

ref 3.  Shoemake, Ken. \"Animating rotation with quaternion curves.\" 
ACM SIGGRAPH computer graphics. Vol. 19. No. 3. ACM, 1985.



Public functions:
qHelp
qCon
qMult
qRot
qqVec
qAngle
qAxis
qtoR
qtoEuler
qDot

-- Package last Modified September 16, 2015, C. Richards

FOR HELP:
Type FunctionName::usage e.g. WaveletSum::usage
Type WaveletHelp[] for sample code
"]
*)


Begin["`Private`"]
(*functions for calculating Euler angles from quaternions*)
SUBqEuler12angles[iR_,iRn_,iRnm_,v3r_]:=(*determine which euler convention is being used and calculate theta 1 and 2 accordingly*)
Module[{th1,th2,group1,group2,group3,group4},
group1=iR=={1,3,1}||iR=={2,1,2}||iR=={3,2,3};
group2=iR=={1,3,2}||iR=={2,1,3}||iR=={3,2,1};
group3=iR=={1,2,1}||iR=={2,3,2}||iR=={3,1,3};
group4=iR=={1,2,3}||iR=={2,3,1}||iR=={3,1,2};
{th1,th2}=If[

group1,
{ArcTan[v3r[[iRn[[1]]]],v3r[[iRnm[[1]]]]],ArcCos[v3r[[iR[[1]]]]]},

If[group2,
{ArcTan[v3r[[iRn[[1]]]],v3r[[iRnm[[1]]]]],-ArcSin[v3r[[iR[[1]]]]]},

If[group3,
{ArcTan[-v3r[[iRnm[[1]]]],v3r[[iRn[[1]]]]],ArcCos[v3r[[iR[[1]]]]]},

If[group4,
{ArcTan[v3r[[iRnm[[1]]]],-v3r[[iRn[[1]]]]],ArcSin[v3r[[iR[[1]]]]]},
{Null,Null}
 ]
]
]
]
]

(*internal function to calculate 3rd euler angle*)
SUBqEuler3rdAngle[Q_,v3r_,theta1_,theta2_,iR_,iRn_]:=Module[{Q1,Q2,Q12,v3n,v3n12,v3nG,magtheta3,m,theta3},
Q1={Cos[theta1/2],0,0,0};
Q1[[1+iR[[1]]]]=Sin[theta1/2];
Q2={Cos[theta2/2],0,0,0};
Q2[[1+iR[[2]]]]=Sin[theta2/2];
Q12=qMult[Q1,Q2];
v3n={0,0,0};(*first zero is added on*)
v3n[[iRn[[3]]]]=1;
v3n12=qRot[Q12,v3n];
v3nG=qRot[Q,v3n];
magtheta3=Abs[ArcCos[v3n12.v3nG]];
m=(Cross[v3n12,v3nG]).v3r;(*to determine sign of angle*)
theta3=Sign[m]*magtheta3
]



(*---------Quaternion algebra functions for handling 3D rotations----------*)
qCon[q_]:=q*{1,-1,-1,-1};

qMult[q1_,q2_]:=Module[{Q1,Q2,product},(*multiply quaternions as vectors. output vector 4*)
Q1=Quaternion[q1[[1]],q1[[2]],q1[[3]],q1[[4]]];
Q2=Quaternion[q2[[1]],q2[[2]],q2[[3]],q2[[4]]];
product=Q1**Q2;
ReplacePart[product,0->List][[;;]]
]

qRot[q_,vec_]:=(*rotates a vector given a quaternion*)Module[{qcon,vec4},
qcon={q[[1]],-q[[2]],-q[[3]],-q[[4]]};
vec4=Join[{0},vec];
qMult[(qMult[q,vec4]),qcon][[2;;4]]
]

qqVec[vec3u_,vec3v_]:=Module[{normUV,cosTheta,halfCos,halfSin,vec3normal},
normUV=Norm[vec3u]*Norm[vec3v];
cosTheta=(vec3u.vec3v)/normUV;
halfCos=Sqrt[0.5*(1+cosTheta)];
halfSin=Sqrt[0.5*(1-cosTheta)];
vec3normal=Cross[vec3u,vec3v]/(normUV*2*halfSin*halfCos);
{halfCos,halfSin*vec3normal[[1]],halfSin*vec3normal[[2]],halfSin*vec3normal[[3]]}
]

qAngle[qUnit_]:=2*ArcCos[qUnit[[1]]];

qAxis[qUnit_]:=Normalize[qUnit[[2;;4]]/Sin[qAngle[qUnit]]]

(*---------Quaternion conversion to Euler angles and Rotation matrices----------*)
qtoR[q_]:=Module[{qx,qy,qz,qw,A,B,CC},
qw=q[[1]];
qx=q[[2]];
qy=q[[3]];
qz=q[[4]];
A={1-2*qy^2-2*qz^2,2*qx*qy-2*qz*qw,2*qx*qz+2*qy*qw};
B={2*qx*qy+2*qz*qw,1-2*qx^2-2*qz^2,2*qy*qz-2*qx*qw};
CC={2*qx*qz-2*qy*qw,2*qy*qz+2*qx*qw,1-2*qx^2-2*qy^2};
{A,B,CC}
]

qtoEuler[Q_,iR_:{3,1,2}]:=Module[{A,iRn,iRnm,v3,v3r,q1,q2,q3},
A=IdentityMatrix[3];(*rotation axes*)
iRn=RotateLeft[iR];
iRnm=RotateLeft[iRn];
v3=A[[iR[[3]]]];
v3r=qRot[Q,v3];
{q1,q2}=SUBqEuler12angles[iR,iRn,iRnm,v3r];
q3=SUBqEuler3rdAngle[Q,v3r,q1,q2,iR,iRn];
{q1,q2,q3}
]

qtoEuler321[q_]:=Module[{qr,qi,qj,qk,eA1,eA2,eA3},
qr=q[[1]];
qi=q[[2]];
qj=q[[3]];
qk=q[[4]];

eA3=ArcTan[1-2*(qi^2+qj^2),2*(qr*qi+qj*qk)];
eA2=ArcSin[2*(qr*qj-qk*qi)];
eA1=ArcTan[1-2*(qj^2+qk^2),2*(qr*qk+qi*qj)];

{eA1,eA2,eA3}

]

(*---------time derivative of quaternion-------------*)
qDot[axisUnitSym_,order_:1]:=Module[{qsym,qdotsym},
qsym={Cos[theta[t]/2],axisUnitSym*Sin[theta[t]/2]};
qdotsym=D[qsym,{t,order}];
Flatten[qdotsym]
]


qHelp[]:=Print["
(*----Example 1.  Plot a frog jumping------*)
xyzPOINTS=sampleDATA[];\[IndentingNewLine]Manipulate[\[IndentingNewLine]videoframe=frame;\[IndentingNewLine]stickFigure=Partition[xyzPOINTS[[frame]],3];(*partition into xyz values*)\[IndentingNewLine]leg=stickFigure[[1;;5]];\[IndentingNewLine]plot1=ListPointPlot3D[leg,PlotRange\[Rule]{{30,130},{-20,80},{-20,80}},BoxRatios\[Rule]{1, 1, 1},AspectRatio\[Rule]1];(*landmark points*)\[IndentingNewLine]plot2=Graphics3D[Line[leg]];(*connect the dots*)\[IndentingNewLine]Show[plot1,plot2],{frame,1,66,1}]

(*----Example 2.  Use a quaternion to rotate a vector------*)
(*we will generate a quaternion to represent the spatial orientation of an arbitrary vector in 3D space with respect to a reference vector*)\[IndentingNewLine](*Step 1. Create a vector and rotate it by 3 arbitrary rotations to create our\[IndentingNewLine]arbitrary vector*)\[IndentingNewLine]vec={0,1,0};(*reference vector*)\[IndentingNewLine]rot1=RotationMatrix[0.5,{1,0,0}];(*rotate by 0.5 radians about x axis*)\[IndentingNewLine]rot2=RotationMatrix[0.8,{0,1,0}];(*rotate by 0.8 radians about y axis*)\[IndentingNewLine]rot3=RotationMatrix[0.2,{0,0,1}];(*rotate by 0.2 radians about y axis*)\[IndentingNewLine]rmatrix=rot1.rot2.rot3;(*composite rotation matrix*)\[IndentingNewLine]vecR=rmatrix.vec;(*rotated vector*)\[IndentingNewLine](*Step 2. represent the rotated vector using a quaternion*)\[IndentingNewLine]q=qqVec[vec,vecR];\[IndentingNewLine]Norm[q](*Verify that it's a unit quaternion*)\[IndentingNewLine](*Step 3.  rotate the reference vector with this quaternion to demonstrate that\[IndentingNewLine]a quaternion can do the same thing as a rotation matrix*)\[IndentingNewLine]vecRalt=qRot[q,vec];\[IndentingNewLine]vecRalt\[Equal]vecR(*verify that the two vectors are the same. In other words, we can use a quaternion to rotate a vector instead of a rotation matrix*)

(*----Example 3.  Extract rotation angle and axis from a unit quaternion------*)
(*make an arbitrary vector whose angle we know by inspection*)\[IndentingNewLine]vec1={1,0,0};(*reference vector*)\[IndentingNewLine]vec2={1,1,0};(*rotated vector*)\[IndentingNewLine]axis=Cross[vec1,vec2](*intuitively we know this is a rotation about the z axis by 45 degrees *)\[IndentingNewLine]q=qqVec[vec1,vec2];\[IndentingNewLine]axisAlt=qAxis[q];\[IndentingNewLine]axisAlt\[Equal]axis(*verify they are the same*)\[IndentingNewLine]angle=qAngle[q](*should be 45 degrees = Pi/4 radians*)

(*-----Example 4.  Convert a quaternion to a rotation matrix*)
(*copy the code from Example 2 then do the following*)
rmatrixAlt=qtoR[q];
rmatrixAlt.vec == rmatrix.vec (*verify you get the same answer*)
(*Note the rotation matrix itself won't necessarily be the same - you can rotate in infinite ways
to arrive at the same final orientation*)

(*-----Example 5.  Convert a quaternion to Euler angles----*)
(*copy the code from Example 2 and do the following*)
eAngles=qtoEuler[q,{1,2,3}];(*get Euler angles for a X-Y-Z rotation*)\[IndentingNewLine](*make the rotation matrices to compose the rotation*)\[IndentingNewLine]rota=RotationMatrix[eAngles[[1]],{1,0,0}];\[IndentingNewLine]rotb=RotationMatrix[eAngles[[2]],{0,1,0}];\[IndentingNewLine]rotc=RotationMatrix[eAngles[[3]],{0,0,1}];\[IndentingNewLine](*pray that you get the same vector*)\[IndentingNewLine]rmatrixAlt=rota.rotb.rotc;\[IndentingNewLine]rmatrixAlt.vec\[Equal]vecR

"]


End[]
EndPackage[]
