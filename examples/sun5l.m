<< "../mbc.m";
<< "../AMBRE.m";

sun5l = Fullintegral[{1}, 
                    {PR[k1, m, n1] 
                     PR[k2, m, n2] 
                     PR[k1 + k2 + k3 + k4 + k5, m, n3] 
                     PR[k3, m, n4] 
                     PR[k4, m, n5] 
                     PR[k5, m, n6]}, 
                    {k5, k4, k3, k2, k1}];
invariants = {};

IntPart[1]
SubLoop[integral]

IntPart[2]
SubLoop[integral]

IntPart[3]
SubLoop[integral]

IntPart[4]
SubLoop[integral];

IntPart[5]
repr = SubLoop[integral];

repr = BarnesLemma[repr, 1];
repr = BarnesLemma[%, 2];

repr = repr /. {n1 -> 1, n2 -> 1, n3 -> 1, n4 -> 1, n5 -> 1, n6 -> 1};

rules = MBoptimizedRules[repr, eps -> 0, {}, {eps}];

integrals = MBcontinue[repr, eps -> 0, rules];

ser = MBexpand[{integrals}, 1/Gamma[1+eps]^5, {eps, 0, -1}];

MBintegrate[ser, {m -> 1}]

(* ******************************************************************
   
   Compare with hep-ph/0403122 eq.41:
   
   3/ep^5 + 16.5/ep^4 + 51.95833/ep^3 + 125.6715/ep^2 + 259.9876/ep +
   347.3551 − 2453.494*ep − 31545.55*ep^2 − 311303.1*ep^3 
   
   ****************************************************************** *)
  
