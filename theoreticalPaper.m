(* ::Package:: *)

BeginPackage["theoreticalPaper`"]

  (*Influence function - Section IV*)
  PlotInfluenceFunctionApprox::usage=" PlotInfluenceFunctionApprox[score_,numC_] plots the influence function for the whole class distribution assuming that the BDR is assesed using 'score' in a 'numC'-class problem (Gaussian mixtures)";
  (*Competitiveness bounds - Section VI*)
  PlotScoreCompetitivenessBounds::usage="PlotScoreCompetitivenessBounds[numC_] plots the competitiveness regions for Holder means used as scores in a 'numC'-class problem";  
  
	Begin["`Private`"]

	PlotInfluenceFunctionApprox[score_,numC_]:=Block[{},
		k = numC;

		(*Bayes Decision Rule Regions; Ri=(a_i,b_i)*)
		fraction[c_,\[Epsilon]_] := If [(1/c-\[Epsilon]/(c-1))== 0, Infinity,(1/c+\[Epsilon])/(1/c-\[Epsilon]/(c-1)) ];

		(* 
			a_i =   -\[Infinity]                                                                        if i=1
            1/(2\[Lambda])*(2*\[Sigma]*Log[(1/k+\[Epsilon])/(1/k-\[Epsilon]/(k-1))]+\[Lambda]^2)                                     if i=2
            max{ \[Lambda]/2*(2i-3) , 1/(2(i-1)\[Lambda])*(2*\[Sigma]*Log[(1/k+\[Epsilon])/(1/k-\[Epsilon]/(k-1))]+(i-1)^2\[Lambda]^2)}    otherwise	 
		*)
		A[\[Epsilon]_?NumericQ,\[Sigma]_,\[Lambda]_?NumericQ,i_]:=Block[{maxVal=+Infinity},
			If[i!= 1,
			If[i==2,maxVal = 1/(2\[Lambda])*(2*\[Sigma]^2*Log[fraction[k,\[Epsilon]]]+ \[Lambda]^2),
				maxVal = Max[1/(2(i-1)\[Lambda])*(2*\[Sigma]^2*Log[fraction[k,\[Epsilon]]]+ (i-1)^2*\[Lambda]^2),(2i-3)*\[Lambda]/2]],
				maxVal=-Infinity];
			maxVal
		];

		(*
			b_i =    \[Infinity]                                        if i=k
            1/(2\[Lambda])*(2*\[Sigma]*Log[(1/k+\[Epsilon])/(1/k-\[Epsilon]/(k-1))]+\[Lambda]^2)     if i=1
            \[Lambda]/2*(2i-3)                                      otherwise	 
         *)
		B[\[Epsilon]_?NumericQ,\[Sigma]_,\[Lambda]_?NumericQ,i_]:=Block[{minVal= -Infinity},
			If[i!= k,
			If[i==1,minVal = 1/(2\[Lambda])*(2*\[Sigma]^2*Log[fraction[k,\[Epsilon]]]+ \[Lambda]^2),
				minVal =(2i-1)*\[Lambda]/2],
				minVal=Infinity];
			minVal
		];

		m::ps="The score `1` is not supported";
		priors[c_,e_] := If[c==1,((1/k)+e),(1/k-e/(k-1))];

		Clear[scoreString];
			scoreString =
			Switch[score,

				"REC",
					"Recall for class 1",
					
				"AMN",
					"Arithmetic Mean among the recalls",

				"ACC",
					"classfication Accuracy",

				"GMN",
					"Geometric Mean among the recalls",
					
				"HMN",
					"Harmonic Mean among the recalls",
				
				"MAX",
					"Maximum Recall",
					
				"MIN",
					"Minimum Recall",
					
				_,
					Message[m::ps,score];Exit[]
			];
		Clear[m];
		
		m[ \[Lambda]_?NumericQ,\[Epsilon]_?NumericQ,k_]:=
			Switch[score,

				"REC",(*Recall_i: (F_i(b_i)-F_i(a_i))*)
					(*Class 1*)i=1;
					Boole[B[\[Epsilon],1,\[Lambda],i]> A[\[Epsilon],1,\[Lambda],i]]*(CDF[NormalDistribution[(i-1)\[Lambda],1],B[\[Epsilon],1,\[Lambda],i]]-CDF[NormalDistribution[(i-1)\[Lambda],1],A[\[Epsilon],1,\[Lambda],i]]),

				"AMN",(*A-mean: (1/K)*Sum_{i=1}^K (F_i(b_i)- F_i(a_i))*)
					(Sum[Boole[B[\[Epsilon],1,\[Lambda],c]> A[\[Epsilon],1,\[Lambda],c]]*(CDF[NormalDistribution[(c-1)\[Lambda],1],B[\[Epsilon],1,\[Lambda],c]]-CDF[NormalDistribution[(c-1)\[Lambda],1],A[\[Epsilon],1,\[Lambda],c]]),{c,1,k}])/k,

				"ACC",(*Accuracy: Sum_{i=1}^K \[Eta]_i(F_i(b_i)- F_i(a_i))*)
					Sum[Boole[B[\[Epsilon],1,\[Lambda],c]> A[\[Epsilon],1,\[Lambda],c]]*priors[c,\[Epsilon]]*(CDF[NormalDistribution[(c-1)\[Lambda],1],B[\[Epsilon],1,\[Lambda],c]]-CDF[NormalDistribution[(c-1)\[Lambda],1],A[\[Epsilon],1,\[Lambda],c]]),{c,1,k}],

				"GMN",(*G-mean: Sqrt[Prod_{i=1}^K (F_i(b_i)- F_i(a_i)),k]*)
					(Product[Boole[B[\[Epsilon],1,\[Lambda],c]> A[\[Epsilon],1,\[Lambda],c]]*(CDF[NormalDistribution[(c-1)\[Lambda],1],B[\[Epsilon],1,\[Lambda],c]]-CDF[NormalDistribution[(c-1)\[Lambda],1],A[\[Epsilon],1,\[Lambda],c]]),{c,1,k}])^(1/k),

				"HMN",(*H-mean: (K)*(Sum_{i=1}^K 1/(F_i(b_i)- F_i(a_i)))^(-1)*)
					If[Product[Boole[B[\[Epsilon],1,\[Lambda],c]> A[\[Epsilon],1,\[Lambda],c]],{c,1,k}]==0,0,k*(Sum[1/(CDF[NormalDistribution[(c-1)\[Lambda],1],B[\[Epsilon],1,\[Lambda],c]]-CDF[NormalDistribution[(c-1)\[Lambda],1],A[\[Epsilon],1,\[Lambda],c]]),{c,1,k}])^(-1)],

				"MAX",(*Max: max[(F_i(b_i)- F_i(a_i))]*)
					Max[Table[Boole[B[\[Epsilon],1,\[Lambda],c]> A[\[Epsilon],1,\[Lambda],c]]*(CDF[NormalDistribution[(c-1)\[Lambda],1],B[\[Epsilon],1,\[Lambda],c]]-CDF[NormalDistribution[(c-1)\[Lambda],1],A[\[Epsilon],1,\[Lambda],c]]),{c,1,k}]],

				"MIN",(*Min: min(F_i(b_i)- F_i(a_i))]*)
					Min[Table[Boole[B[\[Epsilon],1,\[Lambda],c]> A[\[Epsilon],1,\[Lambda],c]]*(CDF[NormalDistribution[(c-1)\[Lambda],1],B[\[Epsilon],1,\[Lambda],c]]-CDF[NormalDistribution[(c-1)\[Lambda],1],A[\[Epsilon],1,\[Lambda],c]]),{c,1,k}]],

				_,
					Message[m::ps,score];Exit[]
			];

		(*Computation of the influence function as expressed in eq. 10 of the manuscript*)
		InfluenceF[\[Epsilon]_?NumericQ,k_]:= NIntegrate[m[\[Lambda],0,k]-m[\[Lambda],\[Epsilon],k],{\[Lambda],0.01,10},AccuracyGoal->3];

		(*Plotting the function for the whole domain of \[Epsilon]*)
		Print[
			Style["Influence Function for the ",18,Black],
			Style[scoreString,18,Black, Bold],Style["assesing the Bayes Decision Rule in a ",18,Black],
			Style[k,18,Black],
			Style["-class problem",18,Black]];
		Print[
			Plot[InfluenceF[\[Epsilon],k],{\[Epsilon],-1/k+0.0001,1-1/k-0.0001},AxesLabel->{"\[Epsilon]","I"}, PlotStyle->Thick,Filling->Axis,
			PlotTheme-> "Detailed",
			PlotRange-> {-4,4},
			PlotLegends-> {"Inf. func."},
			Axes -> True,
			BaseStyle->{FontFamily->"Times",FontSize->20},
			Ticks->{{0,-1/k,0,1-1/k},{-3,0,3}},
			ImageSize->600]
		];
	];

	PlotScoreCompetitivenessBounds[numC_]:=Block[{},
		k = numC;
		(*Bounds, eq. 17 and 18 of the manuscript*)
		SSup[p_]:= ((1/k)*(k-1+1/k^p))^(1/p);
		SInf[p_] := 1/k;

		(*Plotting the function for p in [-50,50]*)
		Print[
			Style["Comparison to the uniformly random classifier in a ",18,Black],
			Style[k,18,Black],
			Style["-class problem",18,Black]];
			p1=
				Plot[{SSup[p],SInf[p]},{p,-50,50},
					Filling->{{(1)  -> 1},{(2)  -> 0}}, 
					AxesLabel->{"p","Score"}, 
					PlotRange->{0,1},
					Epilog->{Style[Text["SUPERIOR",{-25,0.75}],14],
						Style[Text["INFERIOR",{0,0.05}],14],
						Style[Text["INDETERMINATE",{25,0.75}],14]},
					PlotTheme->"Detailed", 
					PlotLegends-> {"Ssup","Sinf"},
					Axes -> True,ImageSize->600
				];
		Print[p1]
	];
	
	End[]

EndPackage[]






