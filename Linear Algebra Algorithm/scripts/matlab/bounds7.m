% Cayley-Menger matrix
% non-symmetric, contains up/low bounds (refined)
%
function B = bounds7()
 
B = [ 0 1 1 1 1 1 1 1;
 1 0 1.161752840 2.337522718 4.64790481 4.6488 4.64704249 2.344450;
 1 1.161752840 0 4.597581148 8.437160981 5.034160367 8.389914898 4.69285569;
 1 2.337522718 4.597581148 0 1.18875409 1.18962649 1.172419875 6.311753550;
 1 4.56206881 7.281876342 1.14554209 0 3.1418 3.14565696 7.080325;
 1 4.56676900 4.899489584 1.174837685 3.07686681 0 3.0837 10.906250;
 1 4.56121449 7.679945620 1.14425809 3.07511296 3.07581444 0 10.652725;
 1 2.319225 4.688165797 6.277208390 5.741850 10.557725 9.977350 0];
%
% end of function bounds7
