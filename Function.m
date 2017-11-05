(*This is a Mathematica binary dump file. It can be loaded with Get.*)                       Windows-x86-64                 ENDCONT Stage 1  S0  Max r Payoff Unequal &CompoundExpression Clear Decrement 
Plus Transpose SetDelayed RandomVariate Greater dT Total 
path RuleDelayed dW stock Pattern "ClearSystemCache cfreal 
Sqrt AppendTo K Do Return Function T If result j &NormalDistribution ValueList 
regr Map Exp current 
Null GreaterEqual 
List All Equal Power R 
Slot Put System` For Global` Set LessEqual Blank Times LaguerreL Increment reg discpayoff  System`Private` First x Fit M callput Simul PatternTest FirstPosition 
Part lastsr Table S n i MapIndexed sigma Stage 2 fN  j"   fF  j   f  f2  À"f:  f\  j   f.  j   f   jz   j&   f  f8  j   j   jD   f,  j0   f  f$  Bfr  ø"f  j   f6   ~H f  ¸"j   f  f  Bj*   f  j    f  jB   f  jR   ft  À"f@  f  f  j|   f  j   fb  f  j   j @ f  j   fT  j>   fx  fP  À&fh  fn f4  jL   fp  fl  BfZ  f^  à"j
   f<  f  ø"fV  j`   f fX Stage 3 ¤bd¢0 B¢ 6¢ : 0 6¢ B
¢ : ¦ R* & (À¤,\hXh&À¤à¤vÈ2& $0"l¤$  J ¢0<$0 0 N *, 0(.¤p0¢0 <¢04$¢p6¢p& P < P (< P J< P v<
¢:<>@B¢08D$J z (BJÄ,*À^ ( ÀZP ¢ÀR¢ÐfJT Æ£ÐHNVX¡|\Ð¤Ð¤È2` ð¤ÜRbj  ¤  (Â ¤,\~h$4¤Df"¤¤0ln¤¤|pÄ,**°L h ($ð¥ü dx (<J x|¢TF<£|~|h ( th$nhÀ¬xð¬`Ä$,*h~Å¤| Ä,*h&Õ¤|`Õä>õìRjõìÜRbjõìür÷íþ>^z¢¤` $ §¤|* °$ (J¦¶¨@vÀ^ JÀÀ
¬ÆÀ0(®ÆÀp°ÆÀ0P²(bJ$th ¤h (À¤,\~ºÀ¤,\~hÂ¤04¾À$,~hÂ¤ 0F¾ÂÂ¤¡¾ÀÄÃ¤¡|¼ÆÃ¬¥r¸È¢ ¤` $h£ ¤|\Ì(JÒ ¦j ($J ("J(lJ(6J (FJ(J§dÖ¶ØÚÜÐÞà°" °0 (ä°|4æ±0ä±|êÀ ,\~ (Á |îÆ $,*~ × ¦|ôÒÆ£Ð,Z (Æ§Ð@¨ø(NJX¦ü¥Z* Æ¤0t2Æ¤µR Æ¬µ@þ ¤   ¬|&&÷íþrÆÀ´V*ÆÀLÆÁ JÆ¶Á0¤Ç¶Á|J "÷ïÿÑ>èìÎðrÊòª¢öúH÷ïÿÑ.â÷ïÿÑ8 8Ô¢Stage 4 d¤(*End of Mathematica binary dump file*) Simul[M_, callput_, n_, T_, K_, R_, sigma_, S0_] := 
    (dT = 1/n; r = R/n; S = Table[0, M, n*T + 1]; S[[All,1]] = S0; 
     dW = Sqrt[dT]*RandomVariate[NormalDistribution[], {M, n*T}]; 
     For[i = 1, i <= n*T, i++, S[[All,i + 1]] = S[[All,i]] + 
        R*dT*S[[All,i]] + sigma*S[[All,i]]*dW[[All,i]]]; Clear[dW]; 
     Payoff[stock_] := If[callput == "Put", Max[K - stock, 0], 
       Max[stock - K, 0]]; cfreal = Table[0, M, n*T]; 
     cfreal[[All,n*T]] = Payoff /@ S[[All,n*T + 1]]; 
     lastsr[path_] := FirstPosition[cfreal[[path]], _?(#1 != 0 & ), {0}][[
       1]]; discpayoff[path_, current_] := If[lastsr[path] != 0, 
       cfreal[[path,lastsr[path]]]*Exp[(-r)*(lastsr[path] - current)], 0]; 
     For[j = n*T - 1, j > 0, j--, reg = {}; Do[If[Payoff[S[[i,j + 1]]] > 0, 
         AppendTo[reg, {S[[i,j + 1]], discpayoff[i, j]}], Null], {i, 1, M}]; 
       regr[x_] = Fit[reg, Exp[-x/2]*{Exp[x/2], LaguerreL[0, x], 
           LaguerreL[1, x], LaguerreL[2, x], LaguerreL[3, x]}, x]; 
       For[i = 1, i <= M, i++, If[Payoff[S[[i,j + 1]]] > 0, 
         If[Payoff[S[[i,j + 1]]] >= regr[S[[i,j + 1]]], cfreal[[i,All]] = 0; 
           cfreal[[i,j]] = Payoff[S[[i,j + 1]]], Null], Null]]]; 
     result = Total[MapIndexed[#1*Exp[(-r)*First[#2]] & , Transpose[cfreal]], 
        2]/M; Clear[cfreal]; ClearSystemCache[]; Return[result])
