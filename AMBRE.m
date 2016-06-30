(*----------------------------------------------------*)
(* AMBRE - Automatic Mellin Barnes REpresentation     *)
(* by K. Kajda                                        *)
(*                                                    *)
(* Published in Comput.Phys.Commun.177:879-893,2007.  *)
(* e-Print: arXiv:0704.2423 [hep-ph]                  *)
(* J. Gluza, K. Kajda, T. Riemann                     *) 
(*----------------------------------------------------*)
(* AMBRE generates Mellin Barnes representations for: *)
(*   - multiloop scalar Feynman integrals             *)
(*   - oneloop n-th order tensor Feynman integrals    *)
(*----------------------------------------------------*)              
Print["by K.Kajda   ver: 1.2.1 \nlast modified 10 Jun 2015\nlast executed on ",
      Apply[Dot,If[StringLength[#]==1,"0"<>#,#]& /@ 
	    Map[ToString[#]&,Reverse[Take[Date[],3]]]],
      " at ",
      Apply[StringJoin,Insert[(If[StringLength[#]==1,"0"<>#,#]& /@ 
	    Map[ToString[#]&,Take[Date[],{4,5}]]),":",2]] ]



BeginPackage["AMBRE`"]


Fauto::usage = "Fauto[0] allows manual modification of F-polynomial. Fauto[1] acitvates automatic mode, set by default"
Fullintegral::usage = "Fullintegral[{numerator},{propagators},{internal momenta}] is a function for input integrals"
invariants::usage = "A list of invariants e.g. {p1^2->s}"
IntPart::usage = "IntPart[iteration] prepares subintegral for given internal momentum."
numerator::usage = "Numerator for a given subloop. Currently it is {1}, no numerators. Please check for future upgrades."
integral::usage = "Propagators for a given subloop"
momentum::usage = "Momentum for a given subloop"
SubLoop::usage = "Subloop[integral]"
BarnesLemma::usage = "BarnesLemma[repr_,1] tries to apply 1st Barnes lemma on a given representation. \n BarnesLemma[repr_,2] is used to apply 2nd Barnes lemma.\n\n These two functions automatically apply lemmas on every M-B integral variable, only if a given variable appears in Gamma functions but not in functions related to kinematcal variables."
PR::usage = "PR[momenta,mass,power of propagator] is a propagator. Mass must be a symbol."
fupc::usage = "F-polynomial variable"
FX::usage = "FX[X[1]+X[2]+...]^2 means (X[1]+X[2]+...)^2"
X::usage = "Feynman parameters of a form X[1],X[2],..."
eps::usage = "Epsilon variable"
z::usage = "M-B integral variable, writen in a form z1,z2,..."
ARint::usage = "M-B integral for tensor diagram. ARint as a function is able to display given integral result: ARint[1,repr] - displays first integral in a result."
g::usage = "Metric tensor"
d::usage = "d - dimension"

Fmanual::mode = "U and F polynomials will be calculated in MANUAL mode. Now you can modify F polynomial (fupc).";
Fauto::mode = "U and F polynomials will be calculated in AUTO mode. In order to use MANUAL mode execute Fauto[0].";


(***********************************************)
(*               Options                       *)
(***********************************************)
Options[IntPart]={Text->True};
Options[SubLoop]={Text->True,Result->False,Xintegration->True};
Options[BarnesLemma]={Text->True,Shifts->False};

Begin["`Private`"]


(***********************************************)
(*               Functions                     *)
(***********************************************)
SFblue[x_]:=StyleForm[x, FontColor -> RGBColor[0, 0, 1]];
SFred[x_]:=StyleForm[x, FontColor -> RGBColor[1, 0, 0]];
ToExprStr[x_,y_]:=ToExpression[ToString[x]<>ToString[y]];
UnsUnApp[x_,y_]:=Reap[Sow[1,Append[x,y]],_,#1&][[2]];
UnsortedUnion[x_]:=Reap[Sow[1, x], _, #1 &][[2]];
ToTake[x_]:=If[Length[x]>iteration,Take[x,iteration],x];
ToList[x_]:=If[Head[x]===Plus,List @@ x,{x}];
ARint[num_,repr_]:=repr[[num]] /.ARint[x_]->x;
IFprint[str__,opt_]:=If[Text/.opt,Print[str]];

OptChange[fun_,new___Rule]:=
  Module[{i},
     Map[
	 If[Length[i=Position[{new},First[#]->_,1,1]]===0, 
	    #,{new}[[Sequence @@ First[i]]]]&, 
	 Options[fun]] ]

(***********************************************)
(*                 Fauto                       *)
(***********************************************)
Fauto[n_]:=
  Block[{},

    Fmode=n;	
    Switch[Fmode,
	   0,
	   Message[Fmanual::mode];
	   fupc=FUpolyfun[integral,momentum,invariants,{Text->True}];
	   Print["fupc = ",fupc],
	   1,
	   Message[Fauto::mode]]
  
 ]


(***********************************************)
(*              Fullintegral                   *)
(***********************************************)
Fullintegral[ARnumerators_,ARintegral_,ARmomentas_]:=
                (FullInt={ARnumerators,ARintegral,ARmomentas});


(***********************************************)
(*             Part of integral                *)
(***********************************************)
IntPart[iter_,options___Rule] :=
  Module[{option,mom,prop,rul1,rul2,(*den,*)a},


	 option=OptChange[IntPart,options]; 

	 mom[j_]:= Last[FullInt][[j]];
	 prop[k_]:= 
	 Times @@ Complement[a[[k]],Flatten[Table[a[[i]],{i,1,k-1}]],{}];

	 Fmode=1;
	 iteration=iter;
         If[iter==1,
	 {
	    znumb={{1,0}};
	    PropList=FullInt[[2]];
	    ToAddList={1}
	 }];
	 
	 rul1 = PR[x_,y_,a_]PR[x_, y_, b_]->PR[x,y,a+b];
	 rul2 = PR[s1__,a_,n1_]PR[s2__,a_,n2_]:>
	         PR[s1,a,n1+n2]/;Expand[(s1)^2-(s2)^2]==0;
	 Propagators=Times @@ DeleteCases[
			Take[PropList,iteration],_Integer] //.{rul1,rul2};

         a = 
	 Table[
	       If[Head[Propagators]===PR,
		   {Propagators},
		   Cases[Propagators,PR[a_.*mom[i]+x_.,__]]],
	     {i,Length[Last[FullInt]]}];
      
    
	 integral=prop[iteration];
	 momentum=mom[iteration];
	 numerator=First[FullInt]/.Times->Dot /.k_^2->k.k;
	 If[numerator[[1]]==1,numerator=1];
	

	 IFprint["numerator=",numerator,
		 "\nintegral=",integral,
		 "\nmomentum=",momentum,option];
       

	 Fauto[Fmode]

 ]


(***********************************************)
(*        Main for given subloop               *)
(***********************************************)
SubLoop[integral_,options___Rule] := 
  Module[{option,subres,sub0,NewProp,NewPropAux,toadd,
	  intover,Subpart,display,result},	 
	 
    
    option=OptChange[SubLoop,options];
   
    IFprint["Iteration nr",iteration,
	    ": >>Integrating over ",momentum,
	    "<<",option];


 (* Used for tadpoles *)

     If[(integral /. PR[__]->PR /.PR^n_.->n)==1 ,
      {
	{mass,lambda} = integral/. PR[_,b_,c_]->{b,c};
	subres = 
	   (-1)^lambda*Gamma[lambda+eps-2]/
	   (Gamma[lambda]*(mass^2)^(lambda+eps-2)); (*dim*)
        znumb=UnsUnApp[znumb,{iteration+1,znumb[[iteration,2]]}];
	Goto[tadpole]
      },{}];

(********************)

     If[Fmode==1,
     {
       IFprint["Computing U & F polynomial in AUTO mode >>Fauto[",Fmode,"]<<",option];
       fupc = Expand[FUpolyfun[integral,momentum,invariants,option]//Simplify]
     },
     {
       IFprint["U & F polynomial was computed by user >>Fauto[",Fmode,"]<<",option]
     }];


	subres = ARintegration[integral, fupc, option];    

Label[tadpole];
 
     PropList = ToTake[PropList];
     ToAddList = ToTake[ToAddList];

     sub0 = subres /.PR[__]->1;

     rul = PR[x_,y_,a_]->PR[x,y,-a];
     NewProp = PowerExpand[subres/sub0]/. PR[x_,y_]^a_->PR[x,y,a]; 
  
       
     If[NewProp==1,{NewPropAux=NewProp*iteration},{},{NewPropAux=NewProp}];


     PropList = UnsUnApp[PropList,NewPropAux /.rul];     
     ToAddList = UnsUnApp[ToAddList,sub0];

     toadd = Times @@ Take[ToAddList,iteration];

     If[iteration==Length[FullInt[[3]]],
     {
       IFprint["Final representation:",option];
       result = toadd*subres;
       If[Head[result]=!=List,
	   {display=SFblue[result]},
	   {If[Result /.option,
	       display=SFblue[result]/. ARint->SFred[ARint],
	       display=SFblue[ARint[SFred[#]]]& /@ 
	         Flatten[Position[result,ARint[x_]]] ]
	   }];
       
       IFprint[display,option]
     },
     {      
       intover = Take[Last[FullInt],iteration];
       intover = #<>##2& @@ (ToString[#]<>"..."& /@ intover);
       IFprint["Representation after integrating over: ",intover,option];
       Subpart = ToExprStr[SubLoop,iteration];
       result = Subpart[SFblue[toadd*sub0],SFred[NewProp]]
     }];

Label["ToResult"];    
       
 Return[result]]


(*************************************************)
(*      Evaluates given integral                 *)
(*************************************************)
ARintegration[integral_,fupc_,option_] :=
  Module[{ipow, Nup,
          tensorial,xint,Xpart,Gpart,
          tempinv,result},


    ipow = List @@ integral /. PR[__,nu_]->nu;
    Nup = Total[ipow]-2+eps; (*dim*)                   

    tensorial = Numk[integral, momentum,numerator,Nup];
    If[Xintegration/.option,
       {xint = xintegral[First[tensorial],integral,ipow];
        Xpart = First[tensorial]*xint /.X[_]->1},
       {Xpart = First[tensorial]*
	        Product[X[i]^(ipow[[i]]-1),{i,1,Length[ipow]}]
       }];

  
    Gpart = Times @@ (Gamma[#]& /@ ipow);

    result = (-1)^(Nup+2-eps)/Gpart*Xpart; (*dim*)

    
    If[Last[tensorial]==0,result = Plus @@ Flatten[result]];   
    If[Last[tensorial]==1,
    {
      tempinv = invariants /.(p_^2->x__)->(p.p->x)
	                   /.(p1_*p2_->x__)->(p1.p2->x);
      
      result = Plus @@ Flatten[result/.tempinv];
    }];

    If[Last[tensorial]>1,
    {
      tempinv=invariants /.(p_^2->x__)->(p.p->x)
	                 /.(p1_*p2_->x__)->(p1.p2->x);
     
      result=Map[ARint[# /.List->Plus]&,result/.tempinv];
    }];

 Return[result]]


(***********************************************)
(*            Mellin-Barnes formula            *)
(***********************************************)
MBarnes[fupc_, integral_, Nup_]:= 
  Module[{imin,fuplist,poly=fupc,nfu,
	  MB1,pow=Nup,MB2,barnes,
	  FXcheck,aux,temp={}},
	   
	
    znumb = Take[znumb,iteration];
    imin = znumb[[iteration,2]];
    
    Label[FXstart];
    fuplist = ToList[Expand[poly]];
 
    nfu = Length[fuplist];
    MB1 = Sum[ToExprStr[z,i+imin],{i,1,nfu-1}]+pow;
    MB2 = Product[fuplist[[i]]^ToExprStr[z,i+imin]*
		  Gamma[-ToExprStr[z,i+imin]],{i,1,nfu-1}];
    barnes = MB2*fuplist[[nfu]]^(-MB1)*Gamma[MB1]/Gamma[pow];
 
    znumb = UnsUnApp[znumb,{iteration+1,imin+nfu-1}];


    FXcheck = Cases[PowerExpand[barnes*aux],FX[_]^_];
    (*Print["barnes = ", barnes];*)
    If[Length[FXcheck] != 0,
       {
	 {poly,pow} = FXcheck /.{FX[x_]^a_}->{x,-a};
	 temp = Append[temp,barnes/. FX[__]->1];
	 imin = znumb[[iteration+1,2]];
	 znumb = Delete[znumb,-1];
         Goto[FXstart]
       }];

    If[Length[temp]!=0,barnes=First[temp]*barnes];


 Return[barnes]]


(***********************************************)
(*            Integration over  x              *)
(***********************************************)
xintegral[fnu_, integral_, ipow_] := 
  Module[{tomult,expr,xonly,numG,denG},
     
      
    tomult = Product[X[i]^ipow[[i]], {i,1,Length[ipow]}];
    expr = PowerExpand[tomult*fnu];
      
    xonly = expr/(expr /.X[_]->1);
    
    numG = xonly /. Power[_,x__]->Gamma[x];
    denG = numG //. Gamma[a_]*Gamma[b_]->Gamma[a+b];
     
       
 Return[numG/denG]]


(***********************************************)
(*         Vector & Tensor parts               *)
(***********************************************)
Numk[integral_, momentum_, numerator_, Nup_] :=
    Module[{num2,rank,FF,P,PQ,AP,
            indices,doexternal,Qin,APobject,
            ext,ind,external,gams,result},
             

    num2=numerator /.momentum.momentum->Sequence[momentum,momentum];
    rank = Length[num2]; 
    
      
    FF[i_]:=MBarnes[fupc,integral,Nup-i];
    PQ[i_]:=Qin /.a_.*p_*X[n_]:>a*p[i]X[n];

    AP[indices_,rank_]:= 
      Module[{poss,rul1,rul2,rul3,rul4,
	      A,result},
    
      	  poss = Select[Permutations[indices],Signature[#]==1&];
       	  rul1 = A[x_]*P[y_]:>Map[A[Take[#,x]]*P[Take[#,-y]]&,poss];
       	  rul2 = P[n_Integer]:>P[ww];
       	  rul3 = A[x_]:>(Times @@ Map[g[Times @@ #]& ,Partition[x,2]]);
       	  rul4 = P[x_]:>(Times @@ Map[P,x]);
    
       	  result =
	      Map[UnsortedUnion[#]&,
	          Table[
		      A[r]P[rank-r],{r,0,rank,2}] 
	              /.rul1 /.rul2 /.rul3 /.rul4];

    Return[result]];

      

    If[rank == 0, 
    {
      result = {Gamma[Nup]*FF[0]};
      Goto[NumkEND];
    }];


(*--- Creating list of indices ---*)
    If[Head[First[num2]]===momentum,
    {
      indices=num2 /.k_[a_]->a;
      doexternal=False
    },
    {
      indices=Table[ToExpression["mu"<>ToString[i]],{i,1,rank}];
      doexternal=True
    }];        
(*--------------------------------*)

        

(*--- Externals onlyif necessary ---*)  
    If[doexternal,
    {
      ext=num2 /.momentum. p_->p;
      ind=Map[II[#]&,indices];
      external=(Apply[Dot,(ext*ind /.p_*II[mu_]->p[mu])]
		 //.momentum[mu_].momentum[nu_]->g[mu*nu])
		  /.Dot->Times
    }];
(*----------------------------------*)
       
      Qin=ToPoly[integral,momentum][[3]];
      APobject=Expand[AP[indices,rank] /.P->PQ]/.Plus->List;
  
      gams=Table[Gamma[Nup-r/2]/(-2)^(r/2)*FF[r/2],{r,0,rank,2}];

       If[doexternal,
	  APobject=Contract[APobject,external]];
     
      result=gams*APobject;
     

     Label[NumkEND];
      
   Return[{result,rank}]]  

(***********************************************)
(*              Contract                       *)
(***********************************************)
Contract[repr_,external_]:=
      Module[{rulA,rulB,rulC,rulD,rulE,result},
    	
    rulA = g[mu1_*mu3_]*g[mu2_*mu3_]:>g[mu1*mu2];
    rulB = g[__]^2->d;
    rulC = g[mu1_*mu2_]*p_[mu2_] -> p[mu1];
    rulD = p1_[x_]*p2_[x_] :> p1.p2/;p1=!=g;
    rulE = p1_[mu1_]^2 :> p1.p1/;p1=!=X;
   
	   
    result = (repr*external) 
	 //.rulA //.rulB //.rulC //.rulD //.rulE;


  Return[result]]


(***********************************************)
(*           Qin, jpol, tt                     *)
(***********************************************)
ToPoly[integral_, momentum_] := 
  Module[{jpol, ttpol, Qin, denom, integrand},

    denom = List @@ integral; Print[Sum[X[i]*denom[[i]],{i,Length[denom]}]];
    integrand =
      Sum[X[i]*denom[[i]] /.PR[k_,m_,_]->k^2-m^2,{i,Length[denom]}];
    ttpol = Coefficient[integrand, momentum*momentum];
    jpol = integrand - momentum*ttpol*momentum;
    Qin = Expand[-1/2*Coefficient[jpol, momentum]];
    jpol = Expand[jpol + 2*momentum*Qin];
    (*Print[{ttpol,jpol,Qin}];*)

 Return[{ttpol,jpol,Qin}]]


(***********************************************)
(*            Finds F function                 *)
(***********************************************)
FUpolyfun[integral_, momentum_, invariants_, option_] := 
  Module[{jpol, ttpol, Qin, (*denom, integrand,*) fpoly, isec},


    {ttpol,jpol,Qin}=ToPoly[integral,momentum];
    Upoly = ttpol;
    fpoly = Expand[Simplify[Expand[Qin*Qin - ttpol*jpol]]]; (*Print[Simplify[ttpol*(Expand[Qin*1/ttpol*Qin] - jpol)]];*)


     If[momentum == Last[FullInt[[3]]],
	 fpoly = Expand[fpoly /. invariants],
	 {},
	 fpoly = Expand[PRfind[fpoly,invariants] /.invariants]
	];

    masslist = Union[List @@ integral /. PR[_, m_, _] -> m];
    masslist = DeleteCases[masslist, 0];
    Do[
  	caslist = Cases[fpoly, masslist[[i]]^2  X[n_]^2];
  	If[Length[caslist] >= 2,
	  {
      	    casexpr = caslist /. a_^2 -> a;
      	    expr = Simplify[Plus @@ casexpr];
      	    FXexpr = expr /. x_.*(y__) -> x^2*FX[y]^2;
      	    FXrule = {First[caslist]->FXexpr-(expr)^2+First[caslist]};
      	    fpoly = Expand[fpoly /. FXrule]; 
      	  },{}],
    {i, 1, Length[masslist]}];

Print["fpoly1 = ", fpoly];

    (***********improvements***************)
    (*negative masses*)
    cleanmassPR = - m_^2 X[i1_] X[i2_] - PR[p_, m_] X[i1_] X[i2_] -> - PR[p, 0] X[i1] X[i2];
    fpoly = fpoly //. cleanmassPR; 
    (*U polynomial*)
    If[Length[fpoly/.{X[_]->1, -1->1}] > 1,{
    (*IF*)
    fxlst = Cases[(List@@fpoly), __ FX[__]^2]; (* Print["fxlst = ", fxlst]; *)
    invlist = DeleteDuplicates[Complement[(List@@fpoly), fxlst] /. {X[_]->1}]; (* Print["invlist = ", invlist]; *)
    
    xlist = Coefficient[Plus@@Complement[(List@@fpoly), fxlst], invlist]; (* Print["xlist = ", xlist]; *)
    
    templist = xlist;
    Do[ If[ Head[templist[[i]]]==Plus,
            templist[[i]] = templist[[i]] /. {Plus->List},{}, 
            templist[[i]] = List[ templist[[i]] ] ],
       {i,1,Length[templist]}]; (* Print["templist1 = ", templist]; *)
    
    Do[ If[ Length[ templist[[i]] ] >= Length[ List @@ Upoly ], 
            If[ Length[ templist[[i]] ] == Length[ List @@ Upoly ],
                {isec = Intersection @@ (templist[[i]] /. {Times->List, X[idx_]^2->List[X[idx], X[idx]]}),
                 If[isec != {}, templist[[i]] = isec]},
                {}
               ] 
           ],
      {i,1,Length[ templist ]}];  (* Print["templist2 = ", templist]; *)

    Do[templist[[i]] = templist[[i]] /. List -> Plus, {i,1,Length[templist]}];
    
    fpoly = Expand[Plus @@ (invlist*templist) + Plus @@ fxlst];
    (*END IF*)
    }];
    (*U polynomial with FX[]*)
(*     If[Length[fpoly/.{X[_]->1, -1->1}] > 1,{
    (*IF*)
    fxlst = Cases[(List@@fpoly), __ FX[__]^2]; Print["fxlst = ", fxlst];
    (* invlist = DeleteDuplicates[Complement[(List@@fpoly), fxlst] /. {X[_]->1}]; Print["invlist = ", invlist]; *)
    invlist = DeleteDuplicates[(List@@fpoly) /. {X[_]->1, FX[_]->1}]; Print["invlist = ", invlist];
    
    (* xlist = Coefficient[Plus@@Complement[(List@@fpoly), fxlst], invlist]; Print["xlist = ", xlist];
    fxlist = Cases[#,FX[__]^2]&/@(List@@##&/@Coefficient[fpoly, invlist]); Print["fxlist = ", fxlist]; *)
    xlist = Coefficient[Plus@@Complement[(List@@fpoly), fxlst], invlist]; Print["xlist = ", xlist];
    fxlist = Coefficient[fpoly, invlist]; Print["fxlist = ", fxlist];

    templist = xlist;
    Do[ If[ Head[templist[[i]]]==Plus,
            templist[[i]] = templist[[i]] /. {Plus->List},{}, 
            templist[[i]] = List[ templist[[i]] ] ],
       {i,1,Length[templist]}]; Print["templist1 = ", templist];
                                      
    Do[ If[ fxlist[[i]] != {}, 
            If[ Length[ templist[[i]] ] == Length[ (fxlist[[i]] /. FX[arg_]^2 -> arg)[[1]] ],
                {tvr = templist[[i]] /. {Times->List, X[idx_]^2->List[X[idx], X[idx]]}, (* Print["tvr = ", tvr], *)
                 isec = Intersection @@ tvr, (* Print["isec = ", isec], *)
                 If[isec != {}, {tocheck = (Plus@@(Complement[#,isec]&/@tvr))[[1]], (* Print["tocheck =", tocheck], *)
                    If[tocheck == (fxlist[[i]] /. FX[arg_]^2 -> arg)[[1]] && tocheck + isec[[1]] == Upoly, 
                       {templist[[i]] = tocheck, fxlist[[i]] = {}}]
                  }]},
                {}
               ] 
           ],
      {i,1,Length[ templist ]}];  

    Do[templist[[i]] = templist[[i]] /. List -> Plus, {i,1,Length[templist]}]; (* Print["templist2 = ", templist]; *)
    Do[fxlist[[i]] = fxlist[[i]] /. {{} -> 0, {var_} -> var}, {i,1,Length[templist]}]; (* Print["fxlist = ", fxlist]; *)

    fpoly = Expand[Plus @@ (invlist*templist) + Plus @@ (invlist*fxlist)];
    (*END IF*)
    }]; *)
    (*linear FX*)
    (*????????*)
    (**************************************)

    IFprint["U polynomial...\n",Upoly,
	    "\nF polynomial...\n",fpoly,option];

 Return[fpoly]

]

(*************************************************)
(*          Creates PRs in F polynomial          *)
(*************************************************)
PRfind[fpoly_,invariants_]:=
  Module[{intmom,masrul,nomass,preprop,p1,p2,ints,intrul,F,
	  PRnom,nom, (*PRmass,*)rul1,rul2,nm,PRadd,aux},
   
     
    intmom = Alternatives @@ (a_. Last[FullInt]);
    masrul = DeleteCases[Union[List @@ integral /.PR[_,a_,_]->a],0];
    nomass =
      Simplify[#]& /@ Collect[fpoly/.(#->0& /@ masrul),X[a_]X[b_]];


    preprop = Plus @@ Cases[nomass+aux,(intmom)^2*_|((intmom)+__)^2*_];
 
	   
    p1 = preprop /.(x_)^2->PR[x];
    p2 = Expand[fpoly-preprop];

    ints= Propagators;
    intrul = ints /.PR[a_,b_,_]->(PR[a]->PR[a,b]+b^2);
    intrul = DeleteCases[List @@ (intrul*(aux->0)),aux->0];

    F = Expand[p1+p2 /.intrul /.invariants];       
    
    PRnom = Cases[F,x_. PR[a_]];
    If[Length[PRnom] == 0,{},
    {
      nom = PRnom /. w_. PR[a_]X[n_]X[m_]->X[n]X[m];
      PRMass = Cases[F,Alternatives @@ (a_. m_^2*nom)];
      rul1 = PR[x_]-m_^2->PR[x,m];
      rul2 =-PR[x_]+m_^2->-PR[x,m];
      nm = Plus @@ Join[PRnom, PRMass];
      PRadd = Collect[nm,nom] /.{rul1,rul2};
      F=F-nm+PRadd
    }];
    
   
 Return[Expand[F /. PR[x_]->PR[x,0]]]]


(*************************************************)
(*          Small Barnes 1 Lemma                 *)
(*************************************************)
BarnesLemma[repres_,num_,options___Rule]:=
    Module[{option, 
	    ZETS, RSC, Zcheck, 
	    BL1, BL2, BL3, BL4, 
	    B2a,
	    Bmes1, Bmes2, Bmes3, (*BLru, pow,*) 
	    Zpowers, ShiftRul, powlist, gamlist,
	    varlist, dimrepr,reprList, BLrulToApply, 
	    before, after,
            counter=0, BList={},zel, denR, ifden, numR, ifnum, result},

    option=OptChange[BarnesLemma,options];
      

    ZETS[zlst_]:=
       Module[{FML,result},
	   FML[x_]:=Flatten[Map[(List @@ #)&, x]];
           result=DeleteCases[FML[FML[zlst]],_Integer]//Union;
       Return[result]];



       RSC[x_,y_]:=ToExprStr[z,#]& /@ 
	  Reverse[Sort[
            ToExpression[StringDrop[ToString[#], 1]]& /@ 
	    Complement[x,y,{}]]];

       Zcheck[x_,z_]:=Cases[x,Gamma[a_.*z+__]|Gamma[a_.*z]];

	  

    (* Using 1st or 2nd Barnes-Lemma *)

       Switch[num,
	      1,
	      
	 BL1[z_]:=
	   Gamma[z+(l1_.)]*Gamma[z+(l2_.)]*Gamma[-z+(l3_.)]*Gamma[-z+(l4_.)]
	   ->Gamma[l1+l3]*Gamma[l1+l4]*Gamma[l2+l3]*
	   (Gamma[l2+l4]/Gamma[l1+l2+l3+l4]);

	 BL2[z_]:=
	   Gamma[z+(l1_.)]^2*Gamma[-z+(l3_.)]*Gamma[-z+(l4_.)] 
	   ->Gamma[l1+l3]^2*Gamma[l1+l4]^2/Gamma[2*l1+l3+l4];

	 BL3[z_]:=
	   Gamma[z+(l1_.)]*Gamma[z+(l2_.)]*Gamma[-z+(l3_.)]^2 
	   ->Gamma[l1+l3]^2*Gamma[l2+l3]^2/Gamma[l1+l2+2*l3];

	 BL4[z_]:=
	   Gamma[z+(l1_.)]^2*Gamma[-z+(l3_.)]^2
	   ->Gamma[l1+l3]^4/Gamma[2*l1+2*l3];

	 Bmes1 = ">> Barnes 1st Lemma will be checked for: ";
	 Bmes2 = ">> Representation after 1st Barnes Lemma:  <<";
	 Bmes3 = "   1st Barnes Lemma was applied for: ";
	 BLrul[z_] := {BL1[z],BL2[z],BL3[z],BL4[z]},
     (*----------------------------------------*)	      
	      2,

	 B2a[z_]:=
	   Gamma[(l1_.)+z]*Gamma[(l2_.)+z]*Gamma[(l3_.)+z]*
	   Gamma[(l4_.)-z]*Gamma[(l5_.)-z]/Gamma[(l6_.)+z]
	   :>Gamma[l1+l4]*Gamma[l2+l4]*Gamma[l3+l4]*Gamma[l1+l5]*
	   Gamma[l2+l5]*Gamma[l3+l5]/(Gamma[l1+l2+l4+l5]*
	   Gamma[l1+l3+l4+l5]*Gamma[l2+l3+l4+l5]) 
	   /;Expand[l1 + l2 + l3 + l4 + l5 - l6] === 0;
 

	 Bmes1 = ">> Barnes 2nd Lemma will be checked for: ";
	 Bmes2 = ">> Representation after 2nd Barnes Lemma:  <<";
	 Bmes3 = "   2nd Barnes Lemma was applied for: ";
	 BLrul[z_] := {B2a[z]}

       ];
      

    (* Preparing representation *)

       reprList=ToList[Expand[repres/.Gamma[x_]:>Gamma[Expand[x]]]];	
    

    (* Creating list of zs NOT to be used *)

       Zpowers=LSTofZ[reprList,1];
      

(*---------------------------*)
(*      Shifting part...     *)
       

       If[(Shifts /. option),
	  ShiftRul=SHIFT[Zpowers];
	  If[Length[ShiftRul]!=0,
	     IFprint[">>",SFred[" Shifting: "],ShiftRul,option];
	     reprList=(#/.ShiftRul/.Gamma[x_]:>
		       Gamma[Expand[x]])& /@ reprList; 
	     Zpowers=LSTofZ[reprList,1],{}]
	  ,{}];

(*---------------------------*)

       powlist=ZETS[Zpowers];


   (* Creating list of zs which appear in Gammas *)

       gamlist=ZETS[LSTofZ[reprList,2]];
	
	       
   (* Creating list of zs to check with B-L*)

       varlist=RSC[gamlist,powlist];   
       dimrepr=Length[gamlist];


   (* Beggining of checking *)

      IFprint[Bmes1,varlist," <<\n   Starting with dim=",dimrepr," representation...",option];
       
       result=reprList;
       Do[
	  zel= varlist[[i]];
	  If[Text/.option,WriteString[$Output,"\n",i,". ","Checking ",zel]];
	  
	  denR = Denominator[First[result]];
	  ifden = Length[Zcheck[denR,zel]];

	  numR = Numerator[First[result]];
	  ifnum = Length[Zcheck[numR,zel]];

	  If[num==1 && ifden!=0,Goto["break"]];
	  If[num==1 && ifnum>4,Goto["break"]];
	  If[num==2 && ifden==0,Goto["break"]];
	  If[num==2 && ifnum>5,Goto["break"]];

	  before=result; BLrulToApply = Dispatch[BLrul[zel]];
	  result=(# /. BLrulToApply)& /@ result; (* Print["BLrulToApply = ", BLrulToApply]; *)
	  after=result;

	  If[before=!=after,
	  {
	    counter=++counter;
	    BList=Prepend[BList,zel];
	    If[Text/.option,WriteString[$Output,"...Barnes Lemma was applied."]]
	  },{}];
	 
  Label["break"],

       {i,1,Length[varlist]}];

      IFprint[Bmes2,option];

       If[Length[BList]==0,
	 IFprint[SFred["   Could not apply Barnes-Lemma"],option],
       {
	 IFprint[SFred[Bmes3],BList,
		 "\n   Obtained representation has: dim=",
		 dimrepr-counter,option]
       }];

 Return[Plus @@ result]]

LSTofZ[xli_,contr_] :=
  Module[{powers,pow,result},
    

     powers=First[FullInt[[2]]]/.PR[_,_,n_]:>Abs[n];
     powers=powers/. Abs[x_]->x;
     powers=#->0& /@ DeleteCases[(List @@ powers),_Integer];


   Switch[contr,
	  1,
     pow=First[xli]/.Gamma[_]->1 /.eps->0/.powers,
          2,
     rul={1/Gamma[x__]^a_.->x,Gamma[x_]^a_.->x};
     pow=Cases[First[xli],Gamma[__]^n_.]/.rul /.eps->0 /.powers
   ];


     result=DeleteCases[
		 If[Head[pow]===Times || Head[pow]===List,
		    List @@ pow,{pow}] /._^x_ :> 
			  Expand[x],_Integer|_Real|_Rational,2];

     If[result==1,result={},{},{}];


 Return[result]]

SHIFT[zlist_] := Module[{SD, prep, crosTab, comel, nocomel, ord, result},
    
    SD[x_] := ToExpression[StringDrop[ToString[x], 1]];
    INT[x_, a_] := 
      Table[If[i == j, {0}, 
          Intersection[Part[x, i], Expand[(a)*Part[x, j]]]], {i, 1, 
          Length[x]}, {j, 1, i}];
    
    	prep = DeleteCases[zlist, _Symbol];

    	crosTab = DeleteCases[INT[prep, -1] // Flatten, 0 | -1, 3];   
    	If[Length[Tab] == 1, GOTO["direct1"], {}];
    
    	comel = DeleteCases[INT[crosTab, 1] // Flatten, 0];
    	If[Length[comel] == 0, GOTO["direct1"], {}];
    
    	nocomel = 
      DeleteCases[
        DeleteCases[#, Alternatives @@ comel] & /@ crosTab, _Symbol];
    	ord = Sort[#, SD[#1] < SD[#2] &] & /@ nocomel;
  
    
    Label["direct1"];
    	result = (First[#] -> First[#] - (Take[#, 1 - Length[#]])) & /@ ord;
    
    Return[result]]

End[]


EndPackage[]
