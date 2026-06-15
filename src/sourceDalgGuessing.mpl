

DalgGuessing := module()

option `Copyright (c) 2026 Bertrand Teguia Tabuguia`, package;

export DalgFunGuess, modDalgFunGuess;

local ordertoktuple, polcheckSol, prevlistnumber, FFixedOrdDegFunGuess, FFixedOrdDegFunGuess2, ADEtoRE, modcheckSol, checkSol,
      deltakdiff, GenMaxlistnumber, MultinomialCount, UnrankMultiset, ADECauchyprod, poch, ADEtermToREterm, modpolcheckSol,FixedOrdDegFunGuess,
      modFixedOrdDegFunGuess, modFFixedOrdDegFunGuess, modFFixedOrdDegFunGuess2,LetGenerateMatrix,LetGenerateIntMatrix;


DalgFunGuess:= proc(L::list,
	      {degADE::posint:=2,
	      degPoly::nonnegint:=2,
       termsOfDegPoly::posint:=1,
	       devars::anyfunc(name):=NULL,
	 startfromord::nonnegint:=0,
	   allPolyDeg::truefalse:=false,
	     sparsity::Or(And(positive, fraction), identical(0)):=0,
	 maxIteration::Or(posint,identical(infinity)):=infinity,
	     approach::identical(recurrence,polynomialsubs):=polynomialsubs,
       inputConstants::set(name):={},
	    linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense,HardSystem):=HardSystem},
		   $)::Or(identical(FAIL),`=`,list(`=`));
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Guessing D-algebraic functions (finding their differential equations)";
		local  Y::anyfunc(name),y::name,x::name,A::anyfunc(name),a::name,n::name,i::nonnegint,
		       c::nothing,N::posint,M::posint,V::list,j::nonnegint,Sinit::list(`=`),k::name,Lf::algebraic,
		       nL::posint:=numelems(L),hasterm::truefalse:=false,ADE::algebraic,RE::algebraic,polEq::algebraic,	
		       Nmax:=ceil(nL/(degPoly+1)),K::list,Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),
		       NDE::algebraic,correct::boolean:=false,Arbconst::list,REcheck::algebraic,NRE::algebraic,NpolEq::algebraic,
		       ADEcheck::algebraic,Param::set:=inputConstants,Meqs,beqs,diffLf,dord;
		#minimal deltak order for starting order startfromord
		N:=binomial(degADE+startfromord,degADE);
		#minimal number of unknown
		M:=(degPoly+1)*N;
		V:=[seq(c[i],i=0..M-1)];
		#ADE variables
		if devars=NULL then 
			`tools/genglobal`('y',{},'reset');
			`tools/genglobal`('x',{},'reset');
			y:=`tools/genglobal`('y');
			x:=`tools/genglobal`('x');
			Y:=y(x)
		else 
			y:=op(0,devars);
			x:=op(1,devars);
			#correctness of the input variables
			if y=x then
				error "invalid input: the variables of the differential equation must be distinct"
			end if;
			Y:=y(x)
		end if;
		#Recurrence variables
		if approach=recurrence then
			a:=`tools/genglobal`('a');
			n:=`tools/genglobal`('n');
			interface(warnlevel=0);
			for j from 0 to degADE-1 do
				cat(k,j):=`tools/genglobal`('k')
			end do;
			K:=[seq(cat(k,j),j=0..degADE-1)];
			interface(warnlevel=4);
			A:=a(n);
			Sinit:=[seq(a(i-1)=L[i],i=1..nL)];
		else
			Lf:=PolynomialTools:-FromCoefficientList(L,x) #add(L[i+1]*x^i,i=0..nL-1)
		end if;
		if M > nL then
			if approach=recurrence then
				if nL-N<termsOfDegPoly*degPoly then
					N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
				end if;
				return ifelse(allPolyDeg,FixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,linsolver,maxIteration,inputConstants),FAIL)
			else
				dord:=PDEtools:-difforder(deltakdiff(Y,x,degADE,N),x);
				diffLf[0]:=Lf;
				for j to dord do 
					diffLf[j]:=diff(diffLf[j-1],x) 
				end do;
				if sparsity<>0 then
					return ifelse(allPolyDeg,FFixedOrdDegFunGuess(Lf,diffLf,dord,degADE,degPoly,Y,N,y,x,linsolver,maxIteration,inputConstants,sparsity),FAIL)
				else
					if nL-N<termsOfDegPoly*degPoly then
						N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
					end if;
					return ifelse(allPolyDeg,FFixedOrdDegFunGuess2(Lf,diffLf,dord,degADE,degPoly,Y,N,y,x,linsolver,maxIteration,inputConstants),FAIL)
				end if
			end if
		end if;
		#initialization - ADE and RE of the first iteration
		ADE:=add(add(V[(degPoly+1)*(j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
		dord:=PDEtools:-difforder(deltakdiff(Y,x,degADE,N),x);
		if approach=recurrence then
			RE:=ADEtoRE(ADE,Y,A,K);
			#write the RE for non-negative indices
			#build the linear system and solve it
			Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=0..M-1)];
			#NegInd: list for substituting terms with negative indices to zero
			NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
			Eq:=subs(NegInd,Eq);
		else	
			diffLf[0]:=Lf;
			for j to dord do 
				diffLf[j]:=diff(diffLf[j-1],x) 
			end do;
			polEq:= normal(eval(ADE,[seq(diff(Y,[x$j])=diffLf[j],j=0..dord)]));
			Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; #[seq(coeff(polEq,x,i),i=0..M-1)] -- not efficient
		end if;
		if linsolver=HardSystem then
			Meqs, beqs := ifelse(approach=recurrence,LetGenerateMatrix(remove(has,Eq,a), V, M),LetGenerateMatrix(Eq, V, M));
			S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
			S:=ifelse(S<>NULL,{seq(V[i] = S[i], i = 1 .. M)},NULL);
		else
			S:=ifelse(approach=recurrence,SolveTools:-Linear(remove(has,Eq,a),V,method=linsolver),
			SolveTools:-Linear(Eq,V,method=linsolver))
		end if;
		if S<>NULL then
			if approach=recurrence then
				if remove(v->rhs(v)=0,S)={} then
					hasterm:=has(Eq,a);
					S:=NULL;
					if hasterm then
						if nL-N<termsOfDegPoly*degPoly then
							N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
						end if;
						return ifelse(allPolyDeg,FixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,linsolver,maxIteration,inputConstants),FAIL)
					end if
				else
					#checking the solution
					REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
				end if
			else
				if remove(v->rhs(v)=0,S)={} then
					S:=NULL
				else
					#checking the solution
					ADEcheck, S, correct:=polcheckSol(S,ADE,diffLf,dord,nL,y,x)
				end if
				
			end if
		end if;
		N:=N+1;
		while not(correct) and N<=Nmax do
			V:=[op(V),seq(c[i],i=M..M+degPoly)];
			NDE:=add(V[M+i]*x^(i-1)*deltakdiff(Y,x,degADE,N),i=1..degPoly+1);
			dord:=PDEtools:-difforder(deltakdiff(Y,x,degADE,N),x);
			ADE:=NDE+ADE;
			if approach=recurrence then
				NRE:=ADEtoRE(NDE,Y,A,K);
				RE:=NRE+RE;
				Eq:=[seq(Eq[i+1]+subs(Sinit,eval(NRE,[n=i,Sum=add])),i=0..M-1)];
				Eq:=[op(Eq),seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=M..M+degPoly+1)];
				NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
				Eq:=subs(NegInd,Eq);
			else
				for j from numelems(diffLf) to dord do 
					diffLf[j]:=diff(diffLf[j-1],x) 
				end do;
				NpolEq:= eval(NDE,[seq(diff(Y,[x$j])=diffLf[j],j=0..dord)]);
				polEq:=polEq+NpolEq;
				Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M+degPoly+1]
			end if;
			M:=M+degPoly+1;
			if linsolver=HardSystem then
				Meqs, beqs := LetGenerateMatrix(Eq, V, M);
				S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
				S:=ifelse(S<>NULL,{seq(V[i] = S[i], i = 1 .. M)},NULL);
			else
				S:=SolveTools:-Linear(Eq,V,method=linsolver)
			end if;
			if S<>NULL then
				if approach=recurrence then
					if remove(v->rhs(v)=0,S)={} or has(map(rhs,S),a) then
						#too few initial values
						hasterm:=has(Eq,a);
						S:=NULL;
						if hasterm then
							if nL-N<termsOfDegPoly*degPoly then
								N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
							end if;
							return ifelse(allPolyDeg,FixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,linsolver,maxIteration,inputConstants),FAIL)
						end if
					else
						#verification step
						REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
					end if
				else
					if remove(v->rhs(v)=0,S)={} then
						S:=NULL
					else
						#checking the solution
						ADEcheck, S, correct:=polcheckSol(S,ADE,diffLf,dord,nL,y,x)
					end if
				
				end if
			end if;
			N:=N+1
		end do;
		if correct then
			ADE:=ifelse(approach=recurrence,subs(S,ADE),ADEcheck);
			if approach=recurrence then
				Param:=Param union {n,op(K)};
				Arbconst:= sort([op(remove(has,indets(REcheck),a) minus Param)])
			else
				Param:=Param union {x,y,seq(diff(Y,[x$i]),i=0..dord)};
				Arbconst:= sort([op(indets(ADEcheck) minus Param)])
			end if;
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				ADE:=subs(Arbconst,ADE)
			end if;
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..dord)},'distributed');
			return ADE=0
		else
			if approach=recurrence then
				if nL-N<termsOfDegPoly*degPoly then
					N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
				end if;
				return ifelse(allPolyDeg,FixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,linsolver,maxIteration,inputConstants),FAIL)
			else 
				if sparsity<>0 then
					return ifelse(allPolyDeg,FFixedOrdDegFunGuess(Lf,diffLf,dord,degADE,degPoly,Y,N,y,x,linsolver,maxIteration,inputConstants,sparsity),FAIL)
				else
					if nL-N<termsOfDegPoly*degPoly then
						N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
					end if;
					return ifelse(allPolyDeg,FFixedOrdDegFunGuess2(Lf,diffLf,dord,degADE,degPoly,Y,N,y,x,linsolver,maxIteration,inputConstants),FAIL)
				end if
			end if
		end if
	end proc:

modDalgFunGuess:= proc(L::list,
	         {degADE::posint:=2,
	         degPoly::nonnegint:=2,
          termsOfDegPoly::posint:=1,
	          devars::anyfunc(name):=NULL,
	    startfromord::nonnegint:=0,
	      allPolyDeg::truefalse:=false,
	        sparsity::Or(And(positive, fraction), identical(0)):=0,
	        approach::identical(recurrence,polynomialsubs):=polynomialsubs,
	    maxIteration::Or(posint,identical(infinity)):=infinity,
	         modulus::posint:=7,
          inputConstants::set(name):={}},
		      $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Guessing D-algebraic functions (finding their differential equations)";
		local  Y::anyfunc(name),y::name,x::name,A::anyfunc(name),a::name,n::name,i::nonnegint,
		       c::nothing,N::posint,M::posint,V::list,j::nonnegint,Sinit::list,k::name,Lf::algebraic,
		       nL::posint:=numelems(L),hasterm::truefalse:=true,ADE::algebraic,RE::algebraic,polEq::algebraic,	
		       Nmax:=ceil(nL/(degPoly+1)),K::list,Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),
		       NDE::algebraic,correct::boolean:=false,REcheck::algebraic,NRE::algebraic,NpolEq::algebraic,diffLf,dord,
		       Terms::list(integer),Meqs,beqs,Aindets,ADEcheck::algebraic,Param::set(name):=inputConstants; #,Arbconst::list
		       
		Terms:= try map(term -> term mod modulus, L) catch : [] end try;
		if Terms = [] or Param <> {} then
			return FAIL
		end if;
		#minimal deltak order for starting order startfromord
		N:=binomial(degADE+startfromord,degADE);
		#minimal number of unknown
		M:=(degPoly+1)*N;
		V:=[seq(c[i],i=0..M-1)];
		#ADE variables
		if devars=NULL then 
			`tools/genglobal`('y',{},'reset');
			`tools/genglobal`('x',{},'reset');
			y:=`tools/genglobal`('y');
			x:=`tools/genglobal`('x');
			Y:=y(x)
		else 
			y:=op(0,devars);
			x:=op(1,devars);
			#correctness of the input variables
			if y=x then
				error "invalid input: the variables of the differential equation must be distinct"
			end if;
			Y:=y(x)
		end if;
		#Recurrence variables
		if approach=recurrence then
			a:=`tools/genglobal`('a');
			n:=`tools/genglobal`('n');
			interface(warnlevel=0);
			for j from 0 to degADE-1 do
				cat(k,j):=`tools/genglobal`('k')
			end do;
			K:=[seq(cat(k,j),j=0..degADE-1)];
			interface(warnlevel=4);
			A:=a(n);
			Sinit:=[seq(a(i-1)=Terms[i],i=1..nL)]
		else
			Lf:=PolynomialTools:-FromCoefficientList(Terms,x) #add(Terms[i+1]*x^i,i=0..nL-1)
		end if;
		#underdetermined system
		if M > nL then
			if approach=recurrence then
				if nL-N<termsOfDegPoly*degPoly then
					N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
				end if;
				return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,maxIteration,maxmodulus,inputConstants),FAIL)
			else
				dord:=PDEtools:-difforder(deltakdiff(Y,x,degADE,N),x);
				diffLf[0]:=modp1(ConvertIn(Lf,x),modulus);
				for j to dord do 
					diffLf[j]:=modp1(Diff(diffLf[j-1]),modulus) 
				end do;
				if sparsity<>0 then
					return ifelse(allPolyDeg,modFFixedOrdDegFunGuess(Lf,diffLf,dord,degADE,degPoly,Y,N,y,x,maxIteration,modulus,inputConstants,sparsity),FAIL)
				else
					if nL-N<termsOfDegPoly*degPoly then
						N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
					end if;
					return ifelse(allPolyDeg,modFFixedOrdDegFunGuess2(Lf,diffLf,dord,degADE,degPoly,Y,N,y,x,maxIteration,modulus,inputConstants),FAIL)
				end if
			end if
		end if;
		#initialization - ADE and RE of the first iteration
		ADE:=add(add(V[(degPoly+1)*(j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
		dord:=PDEtools:-difforder(deltakdiff(Y,x,degADE,N),x);
		if approach=recurrence then
			RE:=ADEtoRE(ADE,Y,A,K);
			#write the RE for non-negative indices
			#build the linear system and solve it
			Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add]) mod modulus) mod modulus,i=0..M-1)];
			#Eq:=map(i -> subs(Sinit, eval(RE, [n = i, Sum = add]) mod modulus) mod modulus, [$0 .. M-1]);
			#NegInd: list for substituting terms with negative indices to zero
			NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
			Eq:=subs(NegInd,Eq);
			Aindets:=[op(indets(Eq,a('integer')))];
			Aindets:=[seq(Aindets[j]=cat(a,j),j=1..numelems(Aindets))];
			Eq:=subs(Aindets,Eq)
		else
			diffLf[0]:=modp1(ConvertIn(Lf,x),modulus);
			for j to dord do 
				diffLf[j]:=modp1(Diff(diffLf[j-1]),modulus) 
			end do;
			polEq:= eval(ADE,[seq(diff(Y,[x$j])=modp1(ConvertOut(diffLf[j],x),modulus),j=0..dord)]) mod modulus;
			#polEq:= eval(ADE,Y=Lf) mod modulus;
			Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M] #[seq(coeff(polEq,x,i),i=0..M-1)]
		end if;
		#solving the linear system
		Meqs, beqs := LetGenerateIntMatrix(Eq,V,M,modulus);
		S:= try [LinearAlgebra:-Modular:-LinearSolve(modulus, <Meqs|beqs>, 1,inplace=false)] catch : NULL end try;
		if S<>NULL then
			S:=convert(S[-1],list); #convert(ifelse(remove(v->v=0,convert(beqs,list))=[],S[-1],S[1]),list);
			if approach=recurrence then
				if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets)) then
					hasterm:=has(Eq,map(rhs,Aindets));
					S:=NULL;
					#too little initial values
					if hasterm then
						if nL-N<termsOfDegPoly*degPoly then
							N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
						end if;
						return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,modulus,inputConstants),FAIL)
					end if
				else
					#checking the solution
					S:=[seq(V[i]=S[i],i=1..M)];
					REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
				end if
			else
				if remove(v->v=0,S)=[] then
					S:=NULL
				else
					S:=[seq(V[i]=S[i],i=1..M)];
					ADEcheck, S, correct:=modpolcheckSol(S,ADE,diffLf,dord,nL,y,x,modulus)
				end if
			end if
		end if;
		N:=N+1;
		while not(correct) and N<Nmax do
			V:=[op(V),seq(c[i],i=M..M+degPoly)];
			NDE:=add(V[M+i]*x^(i-1)*deltakdiff(Y,x,degADE,N),i=1..degPoly+1);
			dord:=PDEtools:-difforder(deltakdiff(Y,x,degADE,N),x);
			ADE:=NDE+ADE;
			if approach=recurrence then
				NRE:=ADEtoRE(NDE,Y,A,K);
				RE:=NRE+RE;
				Eq:=[seq(Eq[i+1]+subs(Sinit,eval(NRE,[n=i,Sum=add]) mod modulus) mod modulus,i=0..M-1)];
				Eq:=[op(Eq),seq(subs(Sinit,eval(RE,[n=i,Sum=add]) mod modulus) mod modulus,i=M..M+degPoly+1)];
				NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
				Eq:=subs(NegInd,Eq);
				Aindets:=[op(indets(Eq,a('integer')))];
				Aindets:=[seq(Aindets[j]=cat(a,j),j=1..numelems(Aindets))];
				Eq:=subs(Aindets,Eq)
			else 
				for j from numelems(diffLf) to dord do 
					diffLf[j]:=modp1(Diff(diffLf[j-1]),modulus) 
				end do;
				NpolEq:= eval(NDE,[seq(diff(Y,[x$j])=modp1(ConvertOut(diffLf[j],x),modulus),j=0..dord)]) mod modulus;
				polEq:=polEq+NpolEq mod modulus;
				Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M+degPoly+1] #[seq(coeff(polEq,x,i),i=0..M+degPoly+1)]
			end if;
			M:=M+degPoly+1;
			Meqs, beqs := LetGenerateIntMatrix(Eq,V,M,modulus);
			S:= try [LinearAlgebra:-Modular:-LinearSolve(modulus, <Meqs|beqs>, 1,inplace=false)] catch : NULL end try;
			if S<>NULL then
				S:=convert(S[-1],list); #convert(ifelse(remove(v->v=0,convert(beqs,list))=[],S[-1],S[1]),list);
				if approach=recurrence then 
					if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets)) then
						hasterm:=has(Eq,map(rhs,Aindets));
						S:=NULL;
						if hasterm then
							if nL-N<termsOfDegPoly*degPoly then
								N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
							end if;
							return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,maxIteration,modulus,inputConstants),FAIL)
						end if
					else
						#verification step
						S:=[seq(V[i]=S[i],i=1..M)];
						REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
					end if
				else
					if remove(v->v=0,S)=[] then
						S:=NULL
					else
						S:=[seq(V[i]=S[i],i=1..M)];
						ADEcheck, S, correct:=modpolcheckSol(S,ADE,diffLf,dord,nL,y,x,modulus)
					end if
				end if
			end if;
			N:=N+1
		end do;
		if correct then 
			ADE:=ifelse(approach=recurrence,subs(S,ADE),ADEcheck);
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..dord)},'distributed');
			return ADE=0
		else
			if approach=recurrence then
				if nL-N<termsOfDegPoly*degPoly then
					N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
				end if;
				return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,maxIteration,modulus,inputConstants),FAIL)
			else
				if sparsity<>0 then
					return ifelse(allPolyDeg,modFFixedOrdDegFunGuess(Lf,diffLf,dord,degADE,degPoly,Y,N,y,x,maxIteration,modulus,inputConstants,sparsity),FAIL)
				else
					if nL-N<termsOfDegPoly*degPoly then
						N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
					end if;
					return ifelse(allPolyDeg,modFFixedOrdDegFunGuess2(Lf,diffLf,dord,degADE,degPoly,Y,N,y,x,maxIteration,modulus,inputConstants),FAIL)
				end if
			end if
		end if
	end proc:


#----------------------------------------------  INTERNAL PROCEDURES ------------------------------------------------

FFixedOrdDegFunGuess2:= proc(  Lf::algebraic,
			   diffLf::table,
			     dord::nonnegint,
			   degADE::posint,
			  degPoly::nonnegint,
				Y::anyfunc(name),
				N::nonnegint,
				y::name,
				x::name,
			linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense,HardSystem),
		     maxIteration::Or(posint,identical(infinity)),
			constants::set(name),
			       $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2026 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,total_perms::nonnegint,
		       nL::posint:=degree(Lf,x)+1,ADE::algebraic,polEq::algebraic,randpick,	
		       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),tl,freqs,val,
		       correct::truefalse:=false,Arbconst::list,ADEcheck::algebraic,ul::list,Meqs,beqs,
		       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint);
		l:=GenMaxlistnumber(N,degPoly,nL-N);
		while l<> FAIL and not(correct) do
			ul:=sort([op({op(l)})]);
			tl:=Statistics:-Tally(l);
			freqs:=Array([seq(eval(val, tl), val = ul)], datatype = integer);
			total_perms:=MultinomialCount(freqs);
			if total_perms < maxIteration then
				Ll:=Iterator:-Permute(l,N);
				M:=add(l)+N;
				for degCoeffs in Ll do
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					polEq:=eval(ADE,[seq(diff(Y,[x$j])=diffLf[j],j=0..dord)]);
					Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; 
					if linsolver=HardSystem then
						Meqs, beqs := LetGenerateMatrix(Eq, V, M);
						S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
						S:=ifelse(S<>NULL,{seq(V[i] = S[i], i = 1 .. M)},NULL)
					else
						S:=SolveTools:-Linear(Eq,V,method=linsolver)
					end if;
					if S<>NULL then
						if remove(v->rhs(v)=0,S)={} then
							S:=NULL
						else
							ADEcheck, S, correct:=polcheckSol(S,ADE,diffLf,dord,nL,y,x)
						end if
					end if;
					if correct then
						break
					end if
				end do
			else
				randpick:=rand(1..total_perms);
				M:=add(l)+N;
				to maxIteration do
					degCoeffs:=UnrankMultiset(randpick(), Array(ul), freqs, N);
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					polEq:=eval(ADE,[seq(diff(Y,[x$j])=diffLf[j],j=0..dord)]);
					Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; 
					Meqs, beqs := LetGenerateMatrix(Eq, V, M);
					S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
					S:=ifelse(S<>NULL,[seq(V[i] = S[i], i = 1 .. numelems(V))],NULL);
					if S<>NULL then
						if remove(v->rhs(v)=0,S)=[] then
							S:=NULL
						else
							ADEcheck, S, correct:=polcheckSol(S,ADE,diffLf,dord,nL,y,x)
						end if
					end if;
					if correct then
						break
					end if
				end do
			end if;
			l:=prevlistnumber(degPoly,l)
		end do;
		if correct then
			ADE:=ADEcheck;
			Arbconst:=sort([op(indets(ADEcheck) 
				minus (constants union {x,y,seq(diff(Y,[x$i]),i=0..dord)}))]);
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				ADE:=subs(Arbconst,ADE)
			end if;
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..dord)},'distributed');
			return ADE=0
		else
			return FAIL
		end if     
	end proc:
	
FFixedOrdDegFunGuess:= proc(   Lf::algebraic,
			   diffLf::table,
			     dord::nonnegint,
			   degADE::posint,
			  degPoly::nonnegint,
				Y::anyfunc(name),
				N::nonnegint,
				y::name,
				x::name,
			linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense,HardSystem),
		     maxIteration::Or(posint,identical(infinity)),
		   inputConstants::set(name),
			 sparsity::fraction,
			       $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2026 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,randpick,
		       nL::posint:=degree(Lf,x)+1,ADE::algebraic,polEq::algebraic,total_combs::posint,	
		       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),pcentge,
		       correct::truefalse:=false,Arbconst::list,ADEcheck::algebraic,nleft,
		       MnL::posint,ZerosV,zV,zzV::list,unkV::list,idx,Meqs,beqs,pool;
		M:=(degPoly+1)*N;
		#the user fixes the maximum of zero-coefficients with minimum possible value M-nL
		pcentge:=max(sparsity,1-nL/M);
		MnL:=ceil(pcentge*M);
		total_combs:=binomial(M,MnL);
		V:=[seq(c[i],i=0..M-1)];
		#userinfo(2,DalgFunGuess,printf("The total combinations are %d \n", total_combs));
		if total_combs < maxIteration then
			ZerosV:=Iterator:-Combination(M, MnL);
			for zV in ZerosV do
				zzV := [seq(c[idx], idx in zV)];
				zzV:=map(t->t=0,zzV);
				unkV:=subs(zzV,V);
				ADE:=add(add(unkV[(degPoly+1)*(j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
				polEq:=eval(ADE,[seq(diff(Y,[x$j])=diffLf[j],j=0..dord)]);
				unkV:=remove(t->t=0,unkV);
				Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..numelems(unkV)]; 
				if linsolver=HardSystem then
					Meqs, beqs := LetGenerateMatrix(Eq, unkV, numelems(unkV));
					S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
					S:=ifelse(S<>NULL,{seq(unkV[i] = S[i], i = 1 .. numelems(unkV))},NULL)
				else
					S:=SolveTools:-Linear(Eq,unkV,method=linsolver)
				end if;
				if S<>NULL then
					if remove(v->rhs(v)=0,S)={} then
						S:=NULL
					else
						ADEcheck, S, correct:=polcheckSol(S,ADE,diffLf,dord,nL,y,x)
					end if
				end if;
				if correct then
					break
				end if
			end do
		else
			to maxIteration do
				pool:=table([seq(i+1=i,i=0..M-1)]);
				nleft:=M;
				if 2*MnL<M then
					zV:=Array(1..MnL);
					for j from 1 to MnL do
						idx:=rand(1..nleft)();
						zV[j]:=pool[idx];
						pool[idx]:=pool[nleft];
						unassign('pool[nleft]');
						nleft:=nleft-1
					end do
				else
					zV:=Array(1..M-MnL);
					for j from 1 to M-MnL do
						idx:=rand(1..nleft)();
						zV[j]:=pool[idx];
						pool[idx]:=pool[nleft];
						unassign('pool[nleft]');
						nleft:=nleft-1
					end do;
					zV:=convert({seq(0..M-1)} minus convert(zV,set),list)
				end if;
				zzV:=[seq(c[idx], idx in zV)];
				zzV:=map(t->t=0,zzV);
				unkV:=subs(zzV,V);
				ADE:=add(add(unkV[(degPoly+1)*(j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
				polEq:=eval(ADE,[seq(diff(Y,[x$j])=diffLf[j],j=0..dord)]);
				unkV:=remove(t->t=0,unkV);
				Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..numelems(unkV)]; 
				Meqs, beqs := LetGenerateMatrix(Eq, unkV, numelems(unkV));
				S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
				S:=ifelse(S<>NULL,{seq(unkV[i] = S[i], i = 1 .. numelems(unkV))},NULL);
				if S<>NULL then
					if remove(v->rhs(v)=0,S)={} then
						S:=NULL
					else
						ADEcheck, S, correct:=polcheckSol(S,ADE,diffLf,dord,nL,y,x)
					end if
				end if;
				if correct then
					break
				end if		
			end do
			
		end if;
		if correct then
			ADE:=ADEcheck;
			Arbconst:=sort([op(indets(ADEcheck) 
				minus (inputConstants union {x,y,seq(diff(Y,[x$i]),i=0..dord)}))]);
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				ADE:=subs(Arbconst,ADE)
			end if;
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..dord)},'distributed');
			return ADE=0
		else
			return FAIL
		end if     
			       
	end proc:


FixedOrdDegFunGuess:= proc(Sinit::list,
			  degADE::posint,
			 degPoly::nonnegint,
			       Y::anyfunc(name),
			       A::anyfunc(name),
			       N::nonnegint,
			       y::name,
			       x::name,
			       a::name,
			       n::name,
			       K::list,
		       linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense,HardSystem),
		    maxIteration::Or(posint,identical(infinity)),
		  inputConstants::set(name),
			      $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
		       nL::posint:=numelems(Sinit),hasterm::truefalse:=false,ADE::algebraic,RE::algebraic,	
		       Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),randpick,val,
		       correct::truefalse:=false,Arbconst::list,REcheck::algebraic,Meqs,beqs,ul::list,tl,
		       l::list(nonnegint),Ll,m::nonnegint,degCoeffs::list(nonnegint),freqs,total_perms;
		
		l:=GenMaxlistnumber(N,degPoly,nL-N);
		ul:=sort([op({op(l)})]);
		tl:=Statistics:-Tally(l);
		freqs:=Array([seq(eval(val, tl), val = ul)], datatype = integer);
		total_perms:=MultinomialCount(freqs);	
		while l<> FAIL and not(correct) do
			if total_perms < maxIteration then
				Ll:=Iterator:-Permute(l,N);
				M:=add(l)+N;
				for degCoeffs in Ll do:
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					RE:=ADEtoRE(ADE,Y,A,K);
					Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=0..M-1)];
					NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
					Eq:=subs(NegInd,Eq);
					if linsolver=HardSystem then
						Meqs, beqs := LetGenerateMatrix(Eq, V, M);
						S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
						S:=ifelse(S<>NULL,{seq(V[i] = S[i], i = 1 .. M)},NULL)
					else
						S:=SolveTools:-Linear(Eq,V,method=linsolver)
					end if;
					if S<>NULL then
						if remove(v->rhs(v)=0,S)={} or has(map(rhs,S),a) then
							hasterm:=has(Eq,a);
							S:=NULL
						else
							REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
						end if
					end if;
					if correct or hasterm then
						break
					end if
				end do
			else
				randpick:=rand(1..total_perms);
				M:=add(l)+N;
				to maxIteration do
					degCoeffs:=UnrankMultiset(randpick(), Array(ul), freqs, N);
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					RE:=ADEtoRE(ADE,Y,A,K);
					Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=0..M-1)];
					NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
					Eq:=subs(NegInd,Eq);
					Meqs, beqs := LetGenerateMatrix(Eq, V, M);
					S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
					S:=ifelse(S<>NULL,[seq(V[i] = S[i], i = 1 .. numelems(V))],NULL);
					if S<>NULL then
						if remove(v->rhs(v)=0,S)=[] or has(map(rhs,S),a) then
							hasterm:=has(Eq,a);
							S:=NULL
						else
							REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
						end if
					end if;
					if correct or hasterm then
						break
					end if
				end do
			end if;
			l:=prevlistnumber(degPoly,l);
			hasterm:=false
		end do;
		if correct then
			ADE:=subs(S,ADE);
			Arbconst:=sort([op(remove(has,indets(REcheck),a) minus ({n,op(K)} union inputConstants))]);
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				ADE:=subs(Arbconst,ADE)
			end if;
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))},'distributed');
			return ADE=0
		else
			return FAIL
		end if     
			       
	end proc:
	
#----------------------------------

modFixedOrdDegFunGuess:= proc(Sinit::list,
			     degADE::posint,
			    degPoly::nonnegint,
				  Y::anyfunc(name),
				  A::anyfunc(name),
				  N::nonnegint,
				  y::name,
				  x::name,
				  a::name,
				  n::name,
				  K::list,
		       maxIteration::Or(posint,identical(infinity)),
			    modulus::posint,
		     inputConstants::set(name),
				 $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
		       nL::posint:=numelems(Sinit),hasterm::truefalse:=false,ADE::algebraic,RE::algebraic,total_perms,	
		       Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),ul,tl,freqs,val,
		       correct::truefalse:=false,REcheck::algebraic,Aindets,Meqs,beqs,randpick,
		       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint); #,Arbconst::list
		
		l:=GenMaxlistnumber(N,degPoly,nL-N);
		ul:=sort([op({op(l)})]);
		tl:=Statistics:-Tally(l);
		freqs:=Array([seq(eval(val, tl), val = ul)], datatype = integer);
		total_perms:=MultinomialCount(freqs);
		while l<> FAIL and not(correct) do
			if total_perms < maxIteration then
				Ll:=Iterator:-Permute(l,N);
				M:=add(l)+N;
				for degCoeffs in Ll do:
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					RE:=ADEtoRE(ADE,Y,A,K);
					Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add]) mod modulus) mod modulus,i=0..M-1)];
					NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
					Eq:=subs(NegInd,Eq);
					Aindets:=[op(indets(Eq,a('integer')))];
					Aindets:=[seq(Aindets[j]=cat(a,j),j=1..numelems(Aindets))];
					Eq:=subs(Aindets,Eq);
					Meqs, beqs := LetGenerateIntMatrix(Eq, V, M, modulus);
					S:= try [LinearAlgebra:-Modular:-LinearSolve(modulus, <Meqs|beqs>, 1,inplace=false)] catch : NULL end try;
					if S<>NULL then
						S:=convert(S[-1],list);
						if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets))  then
							hasterm:=has(Eq,map(rhs,Aindets));
							S:=NULL
						else
							S:=[seq(V[i]=S[i],i=1..M)];
							REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
						end if
					end if;
					if correct or hasterm then
						break
					end if
				end do
			else
				randpick:=rand(1..total_perms);
				M:=add(l)+N;
				to maxIteration do
					degCoeffs:=UnrankMultiset(randpick(), Array(ul), freqs, N);
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					RE:=ADEtoRE(ADE,Y,A,K);
					Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add]) mod modulus) mod modulus,i=0..M-1)];
					NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
					Eq:=subs(NegInd,Eq);
					Aindets:=[op(indets(Eq,a('integer')))];
					Aindets:=[seq(Aindets[j]=cat(a,j),j=1..numelems(Aindets))];
					Eq:=subs(Aindets,Eq);
					Meqs, beqs := LetGenerateIntMatrix(Eq,V,M,modulus);
					S:= try [LinearAlgebra:-Modular:-LinearSolve(modulus, <Meqs|beqs>, 1,inplace=false)] catch : NULL end try;
					if S<>NULL then
						S:=convert(S[-1],list); 
						if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets))  then
							hasterm:=has(Eq,map(rhs,Aindets));
							S:=NULL
						else
							S:=[seq(V[i]=S[i],i=1..M)];
							REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
						end if
					end if;
					if correct or hasterm then
						break
					end if
				end do
			end if;
			l:=prevlistnumber(degPoly,l);
			hasterm:=false
		end do;
		if correct then
			ADE:=subs(S,ADE);
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))},'distributed');
			return ADE=0
		else
			return FAIL
		end if          
	end proc:

modFFixedOrdDegFunGuess2:= proc(Lf::algebraic,
			    diffLf::table,
			      dord::nonnegint,
			    degADE::posint,
			   degPoly::nonnegint,
				 Y::anyfunc(name),
				 N::nonnegint,
				 y::name,
				 x::name,
		      maxIteration::Or(posint,identical(infinity)),
			   modulus::posint,
		    inputConstants::set(name),
				$)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2026 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,total_perms,
		       nL::posint:=degree(Lf,x)+1,ADE::algebraic,polEq::algebraic,val,randpick,	
		       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),ul,tl,freqs,
		       correct::truefalse:=false,ADEcheck::algebraic,Meqs,beqs,
		       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint); #,Arbconst::list
		       
		l:=GenMaxlistnumber(N,degPoly,nL-N);
		while l<> FAIL and not(correct) do
			ul:=sort([op({op(l)})]);
			tl:=Statistics:-Tally(l);
			freqs:=Array([seq(eval(val, tl), val = ul)], datatype = integer);
			total_perms:=MultinomialCount(freqs);
			if total_perms < maxIteration then
				Ll:=Iterator:-Permute(l,N);
				M:=add(l)+N;
				for degCoeffs in Ll do
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					
					polEq:=eval(ADE,[seq(diff(Y,[x$j])=modp1(ConvertOut(diffLf[j],x),modulus),j=0..dord)]) mod modulus;
					Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; 
					Meqs, beqs := LetGenerateIntMatrix(Eq,V,M,modulus);
					S:= try [LinearAlgebra:-Modular:-LinearSolve(modulus, <Meqs|beqs>, 1,inplace=false)] catch : NULL end try;
					if S<>NULL then
						S:=convert(S[-1],list); 
						if remove(v->v=0,S)=[] then
							S:=NULL
						else
							S:=[seq(V[i]=S[i],i=1..M)];
							ADEcheck, S, correct:=modpolcheckSol(S,ADE,diffLf,dord,nL,y,x,modulus)
						end if
					end if;
					if correct then
						break
					end if
				end do
			else
				randpick:=rand(1..total_perms);
				M:=add(l)+N;
				to maxIteration do
					degCoeffs:=UnrankMultiset(randpick(), Array(ul), freqs, N);
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					polEq:=eval(ADE,[seq(diff(Y,[x$j])=modp1(ConvertOut(diffLf[j],x),modulus),j=0..dord)]) mod modulus;	
					Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; #[seq(coeff(polEq,x,i),i=0..M-1)];
					Meqs, beqs := LetGenerateIntMatrix(Eq,V,M,modulus);
					S:= try [LinearAlgebra:-Modular:-LinearSolve(modulus, <Meqs|beqs>, 1,inplace=false)] catch : NULL end try;
					if S<>NULL then
						S:=convert(S[-1],list); 
						if remove(v->v=0,S)=[] then
							S:=NULL
						else
							S:=[seq(V[i]=S[i],i=1..M)];
							ADEcheck, S, correct:=modpolcheckSol(S,ADE,diffLf,dord,nL,y,x,modulus)
						end if
					end if;
					if correct then
						break
					end if
				end do
			end if;
			l:=prevlistnumber(degPoly,l)
		end do;
		if correct then
			ADE:=ADEcheck; 
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..dord)},'distributed');
			return ADE=0
		else
			return FAIL
		end if     
		       
	end proc:
	
modFFixedOrdDegFunGuess:= proc(Lf::algebraic,
			   diffLf::table,
			     dord::nonnegint,
			   degADE::posint,
			  degPoly::nonnegint,
				Y::anyfunc(name),
				N::nonnegint,
				y::name,
				x::name,
		     maxIteration::Or(posint,identical(infinity)),
		          modulus::posint,
		   inputConstants::set(name),
			 sparsity::fraction,
			       $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2026 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
		       nL::posint:=degree(Lf,x)+1,ADE::algebraic,polEq::algebraic,	
		       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),
		       correct::truefalse:=false,ADEcheck::algebraic,Meqs,beqs,
		       MnL::posint,ZerosV,zV,zzV::list,unkV::list,pcentge,total_combs,pool,nleft,idx; #,Arbconst::list
		
		M:=(degPoly+1)*N;
		pcentge:=max(sparsity,1-nL/M);
		MnL:=ceil(pcentge*M);
		total_combs:=binomial(M,MnL);
		V:=[seq(c[i],i=0..M-1)];
		if total_combs < maxIteration then
			ZerosV:=Iterator:-Combination(M, MnL);
			for zV in ZerosV do
				zzV := [seq(c[i], i in zV)];
				zzV:=map(t->t=0,zzV);
				unkV:=subs(zzV,V);
				ADE:=add(add(unkV[(degPoly+1)*(j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
				polEq:=eval(ADE,[seq(diff(Y,[x$j])=modp1(ConvertOut(diffLf[j],x),modulus),j=0..dord)]) mod modulus;
				unkV:=remove(t->t=0,unkV);
				Eq:=PolynomialTools:-CoefficientList(polEq,x);    
				Meqs, beqs := LetGenerateIntMatrix(Eq,unkV,numelems(unkV),modulus);
				S:= try [LinearAlgebra:-Modular:-LinearSolve(modulus, <Meqs|beqs>, 1,inplace=false)] catch : NULL end try;
				if S<>NULL then
					S:=convert(S[-1],list); 
					if remove(v->v=0,S)=[] then
						S:=NULL
					else
						S:=[seq(unkV[i]=S[i],i=1..numelems(S))];
						ADEcheck, S, correct:=modpolcheckSol(S,ADE,diffLf,dord,nL,y,x,modulus)
					end if
				end if;
				if correct then
					break
				end if
			end do
		else
			to maxIteration do
				pool:=table([seq(i+1=i,i=0..M-1)]);
				nleft:=M;
				if 2*MnL<M then
					zV:=Array(1..MnL);
					for j from 1 to MnL do
						idx:=rand(1..nleft)();
						zV[j]:=pool[idx];
						pool[idx]:=pool[nleft];
						unassign('pool[nleft]');
						nleft:=nleft-1
					end do
				else
					zV:=Array(1..M-MnL);
					for j from 1 to M-MnL do
						idx:=rand(1..nleft)();
						zV[j]:=pool[idx];
						pool[idx]:=pool[nleft];
						unassign('pool[nleft]');
						nleft:=nleft-1
					end do;
					zV:=convert({seq(0..M-1)} minus convert(zV,set),list)
				end if;
				zzV:=[seq(c[idx], idx in zV)];
				zzV:=map(t->t=0,zzV);
				unkV:=subs(zzV,V);
				ADE:=add(add(unkV[(degPoly+1)*(j-1)+i+1]*x^i*deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
				polEq:=eval(ADE,[seq(diff(Y,[x$j])=modp1(ConvertOut(diffLf[j],x),modulus),j=0..dord)]) mod modulus;
				unkV:=remove(t->t=0,unkV);
				Eq:=PolynomialTools:-CoefficientList(polEq,x);  
				Meqs, beqs := LetGenerateIntMatrix(Eq,unkV,numelems(unkV),modulus);
				S:= try [LinearAlgebra:-Modular:-LinearSolve(modulus, <Meqs|beqs>, 1,inplace=false)] catch : NULL end try;
				if S<>NULL then
					S:=convert(S[-1],list); 
					if remove(v->v=0,S)=[] then
						S:=NULL
					else
						S:=[seq(unkV[i]=S[i],i=1..numelems(S))];
						ADEcheck, S, correct:=modpolcheckSol(S,ADE,diffLf,dord,nL,y,x,modulus)
					end if
				end if;
				if correct then
					break
				end if
			end do
		end if;
		if correct then
			ADE:=ADEcheck; 
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..dord)},'distributed');
			return ADE=0
		else
			return FAIL
		end if     		       
	end proc:

#-----------------------------------

ordertoktuple := proc(k::posint,n::nonnegint) option remember;
		local m, im, Tkn, j;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		if n=0 then 
		    return [seq(0,j=1..k)]; #Array(1..k, datatype = integer)
		else   
		    Tkn:=ordertoktuple(k,n-1);
		    m:=min(Tkn);
		    im:=ListTools:-Search(m,Tkn);
		    Tkn[im]:=m+1;
		    if im=k then
			return Tkn
		    else
			return [op(Tkn[1..im]),seq(0,j=im+1..k)] #ArrayTools:-Extend(Tkn[1..im],Array(1..k-im,datatype=integer)) 
		    end if     
		end if    
	end proc:
	
deltakdiff := proc(expr,z::name,k::posint:=2,n::nonnegint:=1,$) option remember;
		local tuple;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		tuple:= select(type,ordertoktuple(k,n)-~1,nonnegint);
		return mul(map(d-> diff(expr,[seq(z,1..d)]), tuple))
	end proc:

LetGenerateIntMatrix := proc(Eq::list,V::list,n::integer,m::posint)
		local A, B,i,j;
		A := Matrix(n, n, [ seq([ seq( coeff(Eq[i], V[j]), j=1..n) ], i=1..n) ], datatype=integer);
		B := Vector(n, [ seq( (-subs(map(v -> v=0, V), Eq[i])), i=1..n) ],datatype=integer);
		return LinearAlgebra:-Modular:-Mod(m,A,integer[]),LinearAlgebra:-Modular:-Mod(m,B,integer[])
	end proc:
	
LetGenerateMatrix := proc(Eq::list,V::list,n::integer)
		local A, B,i,j;
		A := Matrix(n, n, [ seq([ seq( coeff(Eq[i], V[j]), j=1..n) ], i=1..n) ]);
		B := Vector(n, [ seq( (-subs(map(v -> v=0, V), Eq[i])), i=1..n) ]);
		return A,B
	end proc:

polcheckSol:= proc(Sol::Or(list,set),
		ADEsol::algebraic,
		diffLf::table,
		  dord::nonnegint,
		    nL::nonnegint,
		     y::name,
		     x::name,
		     $)
	local S::list, ADE::algebraic, j::nonnegint,
	      checkADE::algebraic, deg::extended_numeric;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	S:=map(normal,Sol);
	ADE:=subs(S,ADEsol);
	checkADE:=expand(eval(ADE,[seq(diff(y(x),[x$j])=diffLf[j],j=0..dord)]));
	deg:= ldegree(checkADE,x);
	return ADE, S, evalb(checkADE=0 or deg>=nL-dord)
end proc:

modpolcheckSol:= proc(Sol::Or(list,set),
		   ADEsol::algebraic,
		   diffLf::table,
		     dord::nonnegint,
		       nL::nonnegint,
			y::name,
			x::name,
			m::posint,
			$)
	local ADE::algebraic,j::nonnegint,
	      checkADE::algebraic, deg::extended_numeric;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	ADE:=subs(Sol,ADEsol);
	checkADE:=eval(ADE,[seq(diff(y(x),[x$j])=modp1(ConvertOut(diffLf[j],x),m),j=0..dord)]) mod m;
	checkADE:=Expand(checkADE) mod m;
	deg:= ldegree(checkADE,x);
	return ADE, Sol, evalb(checkADE=0 or deg>=nL-dord)
end proc:

checkSol:= proc(Sol::Or(list,set),
	      REsol::algebraic,
	     NegInd::list,
	      Sinit::list,
		  M::posint,
		 nL::nonnegint,
		  a::name,
		  n::name,
		  $)
	local S::list, RE::algebraic, checkL::list, checkset::set,i::nonnegint;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	S:=map(normal,Sol);
	RE:=subs(S,REsol);
	checkL:=[op(NegInd),op(Sinit)];
	checkset:={seq(normal(subs(checkL,eval(RE,[n=i,Sum=add]))),i=(nL-numelems(Sol)-1)..nL)};
	checkset:=remove(has,checkset,a);
	return RE, S, evalb(checkset in {{0},{}})
end proc:

modcheckSol:= proc(Sol::Or(list,set),
	         REsol::algebraic,
	        NegInd::list,
	         Sinit::list,
		     M::posint,
		    nL::nonnegint,
		     a::name,
		     n::name,
		     m::posint,
		     $)
	local RE::algebraic, checkL::list, checkset::set,i::nonnegint;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	RE:=subs(Sol,REsol);
	checkL:=[op(NegInd),op(Sinit)];
	checkset:={seq(normal(subs(checkL,eval(RE,[n=i,Sum=add]) mod m) mod m),i=(nL-numelems(Sol)-1)..nL)};
	checkset:=remove(has,checkset,a);
	return RE, Sol, evalb(checkset in {{0},{}})
end proc:

prevlistnumber := proc(k::integer, L::list)
	local A, temp_L, n, m, i, j, current_sum, val;

	A := Array(ListTools:-Reverse(L),datatype=integer);
	m:=ArrayNumElems(A);
	n:=0;
	for i from 1 to m do n:=n+A[i] od;

	# Step 1: Find the rightmost index i that we can decrement
	# To keep it a 'partition' (descending), L[i] must remain >= L[i+1]
	for i from m-1 by -1 to 1 do
		if A[i] > 0 and A[i] > A[i+1] then
		    # Decrease this pivot
		    temp_L:=copy(A);
		    temp_L[i] := temp_L[i] - 1;
		    
		    # Step 2: Recalculate the sum used so far
		    current_sum := 0;
		    for j from 1 to i do current_sum:=current_sum+temp_L[j] od;
		    
		    # Step 3: Fill the rest greedily from left to right
		    # We must respect: digit <= k, digit <= previous_digit, and total sum <= n
		    for j from i+1 to m do
			# The digit cannot exceed k, the remaining budget, or the digit to its left
			val := min(k, n - current_sum, temp_L[j-1]);
			temp_L[j] := val;
			current_sum := current_sum + val
		    end do;
		    
		    if current_sum = n then
			return ListTools:-Reverse(convert(temp_L,list))
		    end if
		end if
	end do;

	# If we exit the loop, we've exhausted all partitions for this sum N.
	# Now we must drop the sum to N-1 and start from the Max again.
	if n > 0 then
		return GenMaxlistnumber(m, k, n-1);
	else
		return FAIL # We have reached [0,0,0...]
	end if
end proc:
	
MultinomialCount := proc(counts::Array)
	local total, den, c;
	total := 0;
	for c in counts do total := total + c; od;
	den := 1;
	for c in counts do den := den * factorial(c); od;
	factorial(total) / den;
end proc:
	
UnrankMultiset := proc(R_in, elements::Array, initial_counts::Array, N::integer)
    local R, P, i, j, counts, sub_total;
    R := R_in;
    P := Array(1..N, datatype=integer);
    counts := copy(initial_counts);
    
    for i from 1 to N do
	for j from 1 to ArrayNumElems(elements) do
	    if counts[j] > 0 then
		# If we put elements[j] at position i, 
		# how many unique ways to finish the rest?
		counts[j] := counts[j] - 1;
		sub_total := MultinomialCount(counts);
		
		if R <= sub_total then
		    # This is our element! Move to the next position i.
		    P[i] := elements[j];
		    break
		else
		    # R is further down the list; skip this group
		    R := R - sub_total;
		    # Restore count to try the next unique element for this position
		    counts[j] := counts[j] + 1
		end if
	    end if
	end do
    end do;
    P
end proc:

GenMaxlistnumber := proc(m::posint, k::integer, n::nonnegint)
	local digits, remaining_sum, i, current_digit;

	remaining_sum := n;
	digits := Array(1..m, datatype=integer);

	for i from 1 to m do
		# Greedy choice: take as much as possible, up to k
		current_digit := min(k, remaining_sum);
		digits[i] := current_digit;
		remaining_sum := remaining_sum - current_digit
	end do;

	# Convert list of digits to a single integer
	ListTools:-Reverse(convert(digits,list))
end proc:

ADEtoRE := proc(ADE::algebraic,Y::anyfunc(name),A::anyfunc(name),K::list(name),$)::algebraic;
	local terms, RE;
	option `Copyright (c) 2025 Bertrand Teguia T.`;
	description "convert algebraic differential equations to recurrence equations";
	#Let terms be the terms in the computed differential equation ADE
	if type(ADE,`+`) then
		terms:=[op(ADE)]
	else
		terms:=[ADE]
	end if;
	#Transform each terms of the differential equation to its correspondent for the recurrence equation
	RE:=map(T->ADEtermToREterm(T,Y,A,K),terms);
	#Return the sum recurrence terms
	return add(RE)
end proc:
	
ADEtermToREterm:= proc(term::algebraic,Y::anyfunc(name),A::anyfunc(name),K::list(name),$)::algebraic;
	local x,n,j,q,mterm,Ldiff,Lrec,xpow,c,Cauchyterm,i;
	option `Copyright (c) 2020 Bertrand Teguia T.`;
	description "Conversion of a differential term to a recurrence term.";
	x:=op(Y);
	n:=op(A);
	j:=PDEtools['difforder'](term,x);
	q:=degree(term,diff(Y,[x$j]));
	mterm:=term;
	Ldiff:=[];
	#collect the derivative data iteratively until it remains the monomial: C*x^l
	do
		Ldiff:=[op(Ldiff),`$`(j,q)];
		mterm:=subs(diff(Y,[x$j])^q=1,mterm);
		j:=PDEtools['difforder'](mterm,x);
		q:=degree(mterm,diff(Y,[x$j]))
	until type(mterm,polynom(anything,x));
	xpow:=degree(mterm,x);
	c:=coeff(mterm,x,xpow);
	#Apply the conversion formula
	Lrec:=map(r->poch(n+1,r)*subs(n=n+r,A), Ldiff);
	Cauchyterm:=Lrec[1];
	#Apply the Cauchy product formula for the derivative part
	for i from 2 to numelems(Lrec) do
		Cauchyterm:=ADECauchyprod(Cauchyterm,Lrec[i],n,K[i-1])
	end do;
	return c*subs(n=n-xpow,Cauchyterm)
end proc:

poch := proc(p::algebraic,k::nonnegint,$)::algebraic; 
	local j;
	return mul(p+j,j=0..k-1) 
end proc:

ADECauchyprod := proc(t1::algebraic,t2::algebraic,n::name,k::name,$)::algebraic;
		return Sum(subs(n=k,t1)*subs(n=n-k,t2),k=0..n)
	end proc:
		
end module:

savelib('DalgGuessing',"PathTo/DalgGuessing.mla"):