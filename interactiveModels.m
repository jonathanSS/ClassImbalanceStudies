(* ::Package:: *)

BeginPackage["interactiveModels`"]

  (*l_K related functions - Section VI*)
  Binary::usage="Binary interactive model to visualise the classification error and the geometric mean of the Bayes decision rule applied to a Gaussian mixture model";
  Multiclass::usage="Multi-class interactive model to visualise the classification error and the geometric mean of the Bayes decision rule applied to a Gaussian mixture model";  
  
	Begin["`Private`"]

	(*Binary*)
	Binary := Manipulate[
		(*Parameters of the BINARY model*)
		sigma=1;
		(*Optimal line for the BDR*)
		R:= If[n==0,-Infinity,If[n==1,+Infinity,(2*sigma^2*Log[n/(1-n)]-(0-d)^2+(0+d)^2)/(2((0+d)-(0-d)))]];

		(*Functions to measure the imbalance degree*)
		(*Kullback-Leibler Divergence*)
		KL[n_]:= 
				If[n==0, Return[N[(1-n)Log[2(1-n)]]],
				If[(1-n)==0, Return[N[n Log[2n]]],
				Return[N[n Log[2n]+(1-n)Log[2(1-n)]]]]];
		(*Hellinger Distance*)
		HE[n_]:=N[(1/Sqrt[2])*Sqrt[(Sqrt[n]-Sqrt[1/2])^2+(Sqrt[1-n]-Sqrt[1/2])^2]];
		(*Total Variation Distance*)
		TV[n_]:= N[(1/2)*(Abs[n-1/2]+Abs[1 - n-1/2])];
		(*Chi-Square Distance*)
		CS[n_] := N[2((n-1/2)^2+(1-n-1/2)^2)];
		(*Imbalance Ratio*)
		IR[n_] := If[n/(1-n)<1,(1-n)/n,n/(1-n)];

		Grid[{
		{
			Row[{
				Plot[
					Evaluate@Table[PDF[NormalDistribution[mu,sigma],x],{mu,{0-d,0+d}}],{x,-10,10},
					Epilog->{Thick,Dashed,Line[{{R,-100},{R,100}}]},
					PlotLegends->Placed[{"Class 1","Class 2"},Bottom],
					PlotRange ->All,
					PlotStyle->Thick, 
					Filling-> Axis,
					PlotLabel->"Decision Regions of the BDR",
					ImageSize->{300,200}]
				,
				Plot[
					{KL[n]/KL[1],HE[n]/HE[1],TV[n]/TV[1],CS[n]/CS[1]},{n,0,1},
					Epilog->{Thick,Line[{{n,-100},{n,100}}]},
					PlotLegends->Placed[{"KL","HE","TV","CS"},Bottom],
					PlotRange -> {{0,1},{0,1}},
					PlotStyle->Thick,
					ImageSize->{300,200},
					PlotLabel->"Imbalance Degree"]

				}]
			},
		{
			Row[{
				Plot[
					{1-( n*CDF[NormalDistribution[0-d,sigma],R]+(1-n)*(1 - CDF[NormalDistribution[0+d,sigma],R])),1- Max[n,1-n]},{n,0,1},Epilog->{Thick,Line[{{n,1- (n*CDF[NormalDistribution[0-d,sigma],R]+(1-n)*(1 - CDF[NormalDistribution[0+d,sigma],R]))},{n,1- Max[n,1-n]}}]},PlotRange ->{{0,1},{0,1}},
					PlotLegends->Placed[{"BDR","ZeroR"},Bottom],
					PlotStyle->{{Purple,Thick},{Yellow,Thick}},
					Filling->{1->{2}},
					PlotLabel->"Error of the BDR and ZeroR",
					ImageSize->{300,200}]
				,
				Plot[
					{Sqrt[CDF[NormalDistribution[0-d,sigma],R]*(1 - CDF[NormalDistribution[0+d,sigma],R])],0},{n,0,1},
					Epilog->{Thick,Line[{{n,0},{n,Sqrt[CDF[NormalDistribution[0-d,sigma],R]*(1 - CDF[NormalDistribution[0+d,sigma],R])]}}]},PlotRange ->{{0,1},{0,1}},PlotStyle->{{Purple,Thick},{Yellow,Thick}},
					PlotLegends->Placed[{"BDR","ZeroR"},Bottom],
					Filling->{2->{1}},
					PlotLabel->"G-mean of the BDR and ZeroR",
					ImageSize->{300,200}]
				}]
			},
		{
			TableForm[{
				{Round[KL[n]/KL[1],0.001],
				Round[HE[n]/HE[1],0.001],
				Round[TV[n]/TV[1],0.001],
				Round[CS[n]/CS[1],0.001],
				Round[IR[n],0.001]}},
				TableHeadings-> {
					{Text["Values"]},
					{Text["KL"], Text["HE"], Text["TV"], Text["CS"], Text["IR"]}
				}
			]
		},
		{
			Text["Kullback-Leibler (KL), Hellinger (HE), Total Variation(TV), Chi-Square (CS), Imbalance Ratio (IR)"]}
		}]
		,
		Style["Parameters of the Binary Model",12,Bold],
		{{d,1,"Separation between the means (d)"},0.01,5},
		{{n,0.5,"Prior probability of class 1 (\[Eta])"},0.01,0.99},
		ContentSize->{750,750}]

	(*Multi-class*)
	Multiclass := Manipulate[
		(*Parameters of the MULTICLASS model*)
		Needs["ComputerArithmetic`"];
		(*value not to reach the extrems*)
		t= 0.0000001;
		sigma=1;
		firstmean :=0;
		k= numK;
		means ={};
		Do[AppendTo[means,firstmean + (c-1)*d],{c,1,k}];

		(*Optimal lines for the BDR*)
		fraction[k_,e_] = If [(1/k-e/(k-1))== 0, Infinity,(1/k+e)/(1/k-e/(k-1)) ];

		Clear[A];
		A[e_?NumericQ,sigma_,d_,i_]:=Block[{maxVal=-Infinity, val},
		Do[
			If[a==1, val =( 2 sigma^2 Log[fraction[k,e]]+d^2 ((i-1)^2-(a-1)^2))/(2 (i-a) d),
				val = (d^2 ((i-1)^2-(a-1)^2))/(2 (i-a) d)];
				maxVal=If[val>maxVal,val,maxVal],{a,1,i-1}];
			maxVal
		];
		Clear[B];
		B[e_?NumericQ,sigma_,d_,i_]:=Block[{minVal=Infinity, val},
		Do[
			If[i==1, val =(2 sigma^2 Log[fraction[k,e]]+d^2 ((b-1)^2-(i-1)^2))/(2 (b-i) d),
				val = (d^2 ((b-1)^2-(i-1)^2))/(2 (b-i) d)];
				minVal=If[val<minVal,val,minVal],{b,i+1,k}];
			minVal
		];

		headings ={};
		left = {};
		right ={};
		Do[
			AppendTo[headings,StringJoin["Class ",ToString[c]]];
			AppendTo[left,
				If[A[e0,sigma,d,c]>B[e0,sigma,d,c],"--",Round[A[e0,sigma,d,c],0.001]]];
			AppendTo[right,
				If[A[e0,sigma,d,c]>B[e0,sigma,d,c],"--",Round[B[e0,sigma,d,c],0.001]]];
		,{c,1,k}];

		(*Functions to measure the imbalance degree*)
		(*Kullback-Leibler Divergence*)
		KL[e_,k_]:= Do[If[((k-1)/k-e)>0,Return[((1/k)+e)*Log[1+k*e]+ (((k-1)/k)-e)*Log[1-k*e/(k-1)]],Return[((1/k)+e)*Log[1+k*e]]],{1}];
		(*Hellinger Distance*)
		HE[e_,k_]:= (1/Sqrt[2])*Sqrt[(Sqrt[1/k+e]-Sqrt[1/k])^2+(k-1)(Sqrt[1/k-e/(k-1)]-Sqrt[1/k])^2];
		(*Total Variation Distance*)
		TV[e_,k_]:=Abs[e];
		(*Chi-Square Distance*)
		CS[e_,k_]:=e^2(k+k/(k-1));
		(*Imbalance Ratio*)
		IR[e_,k_] := Max[((1/k)+e)/(1/k-e/(k-1)),(1/k-e/(k-1))/((1/k)+e)];

		priors[c_,e_] := If[c==1,((1/k)+e),(1/k-e/(k-1))];

		Remove[Acc];
		Acc[e1_?NumericQ, sigma_, d_?NumericQ]:=
			Sum[Boole[B[e1,sigma,d,c]> A[e1,sigma,d,c]]*priors[c,e1]*(CDF[NormalDistribution[means[[c]],sigma],B[e1,sigma,d,c]]-CDF[NormalDistribution[means[[c]],sigma],A[e1,sigma,d,c]]),{c,1,k}];

		Remove[Gmean];
		Gmean[e2_?NumericQ, sigma_, d_?NumericQ]:= 
			(Product[Boole[B[e2,sigma,d,c]> A[e2,sigma,d,c]]*
			(CDF[NormalDistribution[means[[c]],sigma],B[e2,sigma,d,c]]-
			CDF[NormalDistribution[means[[c]],sigma],A[e2,sigma,d,c]]),{c,1,k}])^(1/k);

		Grid[{
			{"Decision Regions for the BDR (Ri =(a,b))"},
			{
				TableForm[{left,right},
				TableHeadings-> {
				{Text["a"],Text["b"]},
				headings}]
			},
			{
				Plot[
					Evaluate@Table[priors[mu,e0]*PDF[NormalDistribution[means[[mu]],sigma],x],{mu,1,Length[means]}],{x,-5,40},
					PlotRange ->All,
					PlotStyle->Thick, 
					Filling-> Axis,
					PlotLegends-> means,
					PlotLabel->"Mixture of Gaussians Model",
					ImageSize->{500,400}]
			},
			{
				Row[{
					Plot[
						{1 - Acc[e1, sigma, d], 1-Max[1/k+e1,(1-1/k-e1)/(k-1)]},{e1,-1/k,1-1/k}, PlotRange ->{{-1/k,1-1/k},{0,1}},
						Epilog->{Thick,Line[{{e0,1- Acc[e0,sigma,d]},{e0,1- Max[1/k+e0,(1-1/k-e0)/(k-1)]}}]},
						PlotRange ->{{0,1},{0,1}},
						PlotLegends->Placed[{"BDR","ZeroR"},Bottom],
						PlotStyle->{{Purple,Thick},{Yellow,Thick}},
						Filling->{1->{2}},
						PlotLabel->"Error of the BDR and ZeroR",
						ImageSize->{300,200}]
					,
					Plot[
						{Gmean[e2, sigma, d], 0},{e2,-1/k,1-1/k}, PlotRange ->{{-1/k,1-1/k},{0,1}},
						Epilog->{Thick,Line[{{e0,0},{e0,Gmean[e0, sigma, d]}}]},
						PlotStyle->{{Purple,Thick},{Yellow,Thick}},
						PlotLegends->Placed[{"BDR","ZeroR"},Bottom],
						Filling->{2->{1}},
						PlotLabel->"G-mean of the BDR and ZeroR",
					ImageSize->{300,200}]
				}]
			},
			{ 
				Plot[
						{KL[e,k],HE[e,k],TV[e,k],CS[e,k]},{e,-(1/k),1-1/k},
						PlotLegends->{"KL","HE","TV","CS"},
						Epilog->{Thick,Line[{{e0,0},{e0,100}}]}, 
						PlotRange -> {{-1/k,1-1/k},All},PlotStyle->Thick,ImageSize->{300,400}
				]
			},
		{

			TableForm[{
				{
					Round[KL[e0,k],0.001],
					Round[HE[e0,k],0.001],
					Round[TV[e0,k],0.001],
					Round[CS[e0,k],0.001],
					Round[IR[e0,k],0.001]}
				},
			TableHeadings-> {{Text["Values"]},{Text["KL"],Text["HE"],Text["TV"],Text["CS"],Text["IR"]}}
			]
		},
		{
			Text["Kullback-Leibler (KL), Hellinger (HE), Total Variation(TV), Chi-Square (CS), Imbalance Ratio (IR)"]}
		}]
		,
		Style["\nParameters of the Multi-class Model",12,Bold],
		{{numK,2,"Number of classes (k)"},{2,3,4,5,6,7,8}},
		{{d,1,"Separation between the means (d)"},0.01,5},
		{{e0,0,"Difference among the priors (\[Epsilon])"},-1/k+t,1-1/k-t},
		ContentSize->{750,1200}]

	End[]

EndPackage[]








































