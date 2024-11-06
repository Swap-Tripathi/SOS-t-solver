clc; clear; close; close all;

%% Example 1 (doi: 10.20944/preprints202409.1361.v1): coupled all-to-all connected homogeneous oscillators
%% CAUTION: MIGHT GIVE NEGATIVE SOLUTIONS TO THE FEASIBILITY PROBLEM NEAR THE CURVE IN EXAMPLE 1. CAN BE
%%          MITIGATED BY SOLVING AT SPECIFIC POINTS RATHER THAN THE WHOLE GRID AT ONCE AS DONE IN THIS CODE

Ndatas=10;                      % Plots Ndata*Ndata grid of (alpha,beta) points where stability/synch is checked
aVec=linspace(0,0.5,Ndatas);    % linspace(x1,x2,Ndata) produces Ndata points between x1 and x2 with distance (x2-x1)/(Ndata-1)
bVec=linspace(0,pi/2,Ndatas);   % linspace(x1,x2,Ndata) produces Ndata points between x1 and x2 with distance (x2-x1)/(Ndata-1)
a2Store = [];                   % define a VECTOR (or possibly MATRIX) of variable size
b1Store = [];
global a;                       % global shares a single copy of variable 'a' with all files where "global a" is called
global b;                       
vStore = {};                    % define an ARRAY where elements can be numbers/characters/matrices
wStore = {};
numHarmonics = [];
optimalIterCounter = 1;         % counts number of points where stability/synch is checked. IDEALLY SHOULD START FROM 0.
res = {};                       % result counter: stores when the solution is found or not
for j = 1 : Ndatas
     for i = 1 : Ndatas
          for harmonics = 3 : 2 : 9                 % j:i:k produces a regularly spaced vector j+n*i where n=0,1,2... as long as it is not greater than k

            a = aVec(i);                            % stores the ith point in the "equally spaced points" and assigns it to the global variable a
            b = bVec(j);                            % stores the ith point in the "equally spaced points" and assigns it to the global variable b
            
            dimsX1 = [1,1] * harmonics;             % changes harmonics accoring to line 21
            dimsF  = [3,3];                         % check "solveSDPAlternative.m" for details. 
                                                    % CAN BE AUTOMATED BY CHECKING LARGEST ENTRY OF FIRST AND SECOND COLUMN OF 
                                                    % FF and adding 1 to it
                                                   
            

            [val_four,V,W] = solveSDPAlternative(dimsF, dimsX1);    % checks for lyapunov density and returns infeasible/feasible/solved to val_four
            
            if  strcmp(val_four,'Solved')                       % if val_four is "not" the same as infeasible 
                disp('optimal sol. found')
                a2Store = [a2Store, a];                             % increases the size of the vector by 1 and adds the value of alpha such that feasibility is obtained for (alpha,beta)
                b1Store = [b1Store, b];                             % increases the size of the vector by 1 and adds the value of beta such that feasibility is obtained for (alpha,beta)
                numHarmonics = [numHarmonics, harmonics];           % stores the harmonic in the line 21 for which solution is obtained
                vStore{optimalIterCounter} = V;                     % if the solution exists, stores V to the vstore ARRAY
                wStore{optimalIterCounter} = W;                     % if the solution exists, stores W to the vstore ARRAY
                res{optimalIterCounter} = val_four;                 % stores feasible/solved to result ARRAY
                optimalIterCounter = optimalIterCounter +1;         % increments optimal counter by 1 any time i or j (alpha or beta) is changed
                break;
            else
                disp('infeasible')
            
            end
    
          end
     end
end

xxxx=linspace(0,pi/2,Ndatas);                                       % the beta vector 
yyyy=0.5*cos(xxxx);                                                 % corresponding vector of evaluations (WHY COS?) ans: this is the known curve for 3-dimensional systems
plot(b1Store,a2Store,'*',xxxx,yyyy)                                 % plots "*" for feasible points (alpha,beta) and line plot yyyy w.r.t xxxx


