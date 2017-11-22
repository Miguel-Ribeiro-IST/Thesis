(*This is a Mathematica binary dump file. It can be loaded with Get.*)         
             Windows-x86-64                 ENDCONT Stage 1 �Equal System` Simul Global` Times Payoff reg dT stock For j 
Slot FirstPosition GreaterEqual &NormalDistribution ValueList LaguerreL Blank K MapIndexed  LessEqual result PatternTest Set T Increment dW Power &CompoundExpression Put S0 
path Table Unequal cfreal Max All Exp callput r  System`Private` i RandomVariate 
regr R 
Sqrt RuleDelayed x Transpose Clear Greater Return discpayoff Decrement First Pattern SetDelayed Map lastsr Total current 
Part If Function sigma 
Plus AppendTo Do 
Null 
List Fit "ClearSystemCache M n S Stage 2 �b  �B   r  �N  �&~  �x  �B:  �"&   @   V  �d   �   Z  �<  �  � ��    n   8   T   � �L  ��  �  �z ��   \   0  �j  ��    �� 6  �"  �&H   2  �B$  ��  ��  ��  �   p  �   
  �"�   �  �X   f  �J  �"�   � �`  �"D �.   �  �R   4   F  �l  �,  �h  �  ��  �( �^   |   �  �" �B   v  �  �  �Stage 3 ����� t �����V��������B *����B P����B t����B �
����d
��$Vp����� ������ P����.4\Ф��T����6𬈁�����n> ���H����X"����n"����"����p"����"�����"����b"����"��� $&(*,.02���� �����V P6����F"8����V�6����F(<����� P����h *X@����F4B����.4, P����FF�$f" ������dX����LN����VJP����F&R����F\ P��$v\�@\�����\ P����.4,\���.4,\¤�V�"`�$.&,\¤�Vb`d¤��`bfä��F^hì���VXZj�x&����>"���Tp����V �>�����t����`v *����V ������>z����`| *����l x~����r�����h *X����FD�Ơ$.D, ��.4,@Ҡ��2T�נ��F�������"�����Ā.D�������0 P����r� *����L�����8"�����d *Ƅ��z���Ƅ��.� Pƌ��������<"���$����������r� *Ā.D������V �(����V �<�����������V�������Ƥ��V��Ƥ��l�� *Ƭ��������F��Ȁ: *���R����FP�Фȁ: *��$\����d�ई�|P��܁l��J����d\ PX���*������"���6�����V �����V����������V�����NP�����F���� v\X�$.D\,Ť��F� *Ą.D\դ��F��䈁���윁l�J��܁l��J�����V�Z������������������������0 t�����ƀ��V��Ƅ���Ƅ��V��Ƅ��L��^DƆ��~��Ɔ��
� t����X �Ɔ��V��ǆ��Fj��xD ��Z$tj"����:>DHTln����������������4��������Stage 4  �(*End of Mathematica binary dump file*) Simul[M_, callput_, n_, T_, K_, R_, sigma_, S0_] := 
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
 
M = 100
 
callput = "Put"
 
n = 12.
 
T = 1.
 
i = 9
 
x = 4
