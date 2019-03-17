% Top Level Script which sets some parameters and calls each method
% After all the methods are run, performance function perf is called to 
% draw a comparison plot

n=10;% dimension of problem
N=10; % number of iterations
x0=10*rand(n,1)-5;% start points same for each method

fprintf('\nNewton\n');
Newton
ElapsedTime = tElapsed';
Iterations = counter';

fprintf('\nBFGS\n');
BFGS
ElapsedTime = horzcat(ElapsedTime,tElapsed');
Iterations = horzcat(Iterations,counter');

fprintf('\nInverse BFGS\n');
InverseBFGS
ElapsedTime = horzcat(ElapsedTime,tElapsed');
Iterations = horzcat(Iterations,counter');

fprintf('\nDFP\n');
DFP
ElapsedTime = horzcat(ElapsedTime,tElapsed');
Iterations = horzcat(Iterations,counter');

fprintf('\nInverse DFP\n');
InverseDFP
ElapsedTime = horzcat(ElapsedTime,tElapsed');
Iterations = horzcat(Iterations,counter');

% fprintf('\nSteepest Descent\n');
% SteepestDescent
% ElapsedTime = horzcat(ElapsedTime,tElapsed');
% Iterations = horzcat(Iterations,counter');

fprintf('\nAll Methods run. DONE!!!!\n');

perf(ElapsedTime,N);