<< "../mbc.m";
<< "../AMBRE.m";

tad3l = Fullintegral[{1}, 
                    {
                        PR[k1, m, 1] 
                        PR[k2, m, 1]
                        PR[k3, m, 1] 
                        PR[k1+k2+k3, m, 1] 
                    }, 
                    {k1, k2, k3}];
invariants = {};

IntPart[1]
SubLoop[integral]

IntPart[2]
SubLoop[integral]

IntPart[3]
repr = SubLoop[integral]

rules = MBoptimizedRules[repr, eps -> 0, {}, {eps}];

integrals = MBcontinue[repr, eps -> 0, rules];

ser = MBexpand[{integrals}, 1/Gamma[1+eps]^3, {eps, 0, 1}];

resmb = MBintegrate[ser, {m -> 1}, MaxNIntegrateDim->1];

Print["Result = ",resmb];

(* ******************************************************************
   
   Compare with hep-ph/0102033 eq.198:
   
   2/ep^3 + 7.666666666667/ep^2 + 17.5/ep +
   22.91666666667 + 21.25179105129*ep − 184.2300051053*ep^2 − 
   661.1105861534*ep^3 − 3685.054779382*ep^4
   
   ****************************************************************** *)
  
