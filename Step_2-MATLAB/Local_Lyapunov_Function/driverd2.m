clc; clear; close; close all;


%% example for N  = 3
aVecT=[];
bVecT=[];
Ndata=10;
aVec=linspace(0,0.5,Ndata);
bVec=linspace(0,pi/2,Ndata);
a2Store = [];
b1Store = [];
r=1; % To find stability in the region Diameter<pi/2r
global a;
global b;
elapsedTime = 0;
totalTime = 0;
averageElapsedTime = 0;
VStore = {};
WStore = {};
XStore = {};
YStore = {};

optimalIterCounter = 1;
res = {};

for j = 1 : Ndata

  for i = 1 : Ndata

    disp(['Estimated Time Remaining: ', num2str(averageElapsedTime * (Ndata*Ndata - ((j-1) * Ndata + i))),' seconds'])

    tic;
    a = aVec(i);
    b = bVec(j);  


    dimsGSV = [5, 5];
    dimsF   = [3, 3];
    
    

    [val_four, V, W, X, Y] =...
        solveSDPAlternative(dimsF, dimsGSV, r);

   
    if strcmp(val_four,'Solved') 
        disp('optimal sol. found')
        a2Store = [a2Store, a];
        b1Store = [b1Store, b];

        VStore{optimalIterCounter}  = V;
        WStore{optimalIterCounter} = W;
        XStore{optimalIterCounter} = X;
        YStore{optimalIterCounter} = Y;

       
        res{optimalIterCounter} = val_four;
        optimalIterCounter = optimalIterCounter +1;
    else
        disp('infeasible')
    end


        elapsedTime = toc;
        totalTime = totalTime + elapsedTime;
        averageElapsedTime = totalTime / ((j-1) * Ndata + i);
    
  end

end

xxxx=linspace(0,pi/2,10);
yyyy=0.5*cos(xxxx);
plot(b1Store,a2Store,'*',xxxx,yyyy)


