tElapsed = zeros(1,N); % Time taken for each outer loop
counter = zeros(1,N);% Inner Loop counter for each outer iteration
for i = 1:N
    tStart = tic;  % Start Time counter for i
    fprintf('Iteration Number = %i\n',i);
    x=x0;
    [f,g] = obj(x);% Evaluating values of function and gradient
    
    k = 0;
    alpha=[];
    K=abs(f);
    tol=sqrt(eps);
    df=1;
    ndx=1;
    
    count = 0;
    while norm(g)>tol && df>100*K*eps && ndx>tol
        count = count + 1;
        p = -g; % Steepest Descent
        [alpha]=ls_V2(k,x,p,alpha);% Line Search algorithm
        dx = alpha*p;
        ndx = norm(dx);
        x = x+dx;
        
        [fnew,gnew] = obj(x);
        df = abs(fnew-f);
        
        p = -gnew;
        if p'*gnew>=0
            fprintf('Not a descent direction, p''*g: %g\n', p'*gnew)
            disp('Finished processing this problem.');
            break
        end
        f = fnew;
        g = gnew;
        k = k+1;
    end
    tElapsed(i) = toc(tStart);  % Stop Time counter
    counter(i) = count;
end