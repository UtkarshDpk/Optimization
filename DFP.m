tElapsed = zeros(1,N); % Time taken for each outer loop
counter = zeros(1,N);% Inner Loop counter for each outer iteration
for i = 1:N
    tStart = tic;  % Start Time counter for i
    fprintf('Iteration Number = %i\n',i);
    x=x0;
    [f,g] = obj(x);% Evaluating values of function & gradient
    
    k = 0;
    alpha=[];
    K=abs(f);
    tol=sqrt(eps);
    df=1;
    ndx=1;
    
    h=eye(n);% Initialising hessian h as Identity matrix
    count = 0;
    while norm(g)>tol && df>100*K*eps && ndx>tol
        count = count + 1;
        p = -h\g;% Descent Direction for DFP
        [alpha]=ls_V2(k,x,p,alpha);% Line Search algorithm
        dx = alpha*p;
        ndx = norm(dx);
        x = x+dx;
        
        [fnew,gnew] = obj(x);
        df = abs(fnew-f);
        
        Y = (gnew -g);
        if(dx'*Y>=0)
            % Calculating new h for DFP
            gama = 1/(dx'*Y);
            % gama is a scalar
            A = eye(n) - gama*(dx*Y');
            hnew = A'*h*A + gama*(Y*Y');
            p = -hnew\gnew;% New Descent Direction
            if (p'*gnew>=0)
                fprintf('Not a descent direction, p''*g: %g\n', p'*gnew);
                fprintf('Finished processing this iteration\n');
                break
            end
        else
            fprintf('Condition dx''*Y>0 not satisfied dx''*Y: %g\n',dx'*Y);
        end
        f = fnew;
        g = gnew;
        h = hnew;
        k = k+1;
    end
    tElapsed(i) = toc(tStart);  % Stop Time Counter
    counter(i) = count;
end