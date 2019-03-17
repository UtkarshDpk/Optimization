function [alpha,counter]=ls_V2(k,x,p,alpha)
% LS_MAT
% finds an alpha value satisfying the Wolfe conditions
% Input parameters: the current point x and the search direction p.
% Author J Kinsella
% Date 27 January 2014
% The code is based on the pseudo-code in the paper
% Algorithm 851: CG DESCENT, a Conjugate Gradient Method
% with Guaranteed Descent
% by   William W. Hager and Hongchao Zhang
% in ACM Transactions on Mathematical Software,
% Vol. 32, No. 1, March 2006, Pages 113--137.
% Available at http://clas.ufl.edu/users/hager/papers/CG/cg_compare.pdf
% Corrected edition 18 March 2014
% Final (?) version 24 March 2014


verb=0; %no debugging output

c1=0.01;
c2=0.4;
rho=5;%1.5;
gamma=0.66;
QuadStep=true;
epsk=sqrt(eps); % temp setting - updated by cgm code
c=init(k,x,QuadStep,alpha);
[a,b]=bracket(@phi,c,rho,epsk);
done=0;
counter=0;
while ~done
    if abs(b-a)<=100*eps
        done=true;
    else
        [anew,bnew]=secant2(@phi,a,b,epsk);
        if (bnew-anew)>gamma*(b-a)
            c=(anew+bnew)/2;
            [anew,bnew]=cg_update(@phi,anew,bnew,c,epsk);
        end
    end;
if ~done
    a=anew;
    b=bnew;
end
    W= wolfe(@phi, c1,c2,c,epsk);
    done=W|| done;%abs(a-b)<=100*eps;
    counter=counter+1;
end
alpha=c; % why not alpha=c??
if ~wolfe(@phi, c1,c2,c,epsk)
    fprintf('FAILED.\n')
    error('failed wcs')
end



% Begin search for alpha sat SW2


[phi0,dphi0]=phi(0,x,p);

[phialpha,dphialpha]=phi(alpha,x,p);
[W,W1,W2,Wapp,SW2]=WOLFEfast(alpha,c1,c2,phi0,dphi0,phialpha,dphialpha,epsk);
done=W&&SW2;
sw2count=0;
sw2max=10;
A=0;
B=alpha;
phiA=phi0;
dphiA=dphi0;
phiB=phialpha;
dphiB=dphialpha;

while ~done && sw2count<=sw2max
    if verb, fprintf('SW2 fails@ alpha: %g.\n',alpha),end
    if verb, fprintf('counter: %i.\n',sw2count),end
    [r1,r2,phir1,dphir1,phir2,dphir2]=invcubrat(@phi,A,B,x,p,phiA,dphiA,phiB,dphiB);
    
    if r1<0 && r2<0
        error('inv cubic interp failed')
    end
    
    if r1>=0
        [W,W1,W2,Wapp,SW2]=WOLFEfast(r1,c1,c2,phi0,dphi0,phir1,dphir1,epsk);
        if W&&SW2
            alpha=r1;
            if verb,      fprintf('r1 used: %g.\n',r1),end
            done=1;
        end;
    end;
    if r2>=0
        [W,W1,W2,Wapp,SW2]=WOLFEfast(r2,c1,c2,phi0,dphi0,phir2,dphir2,epsk);
        if W&&SW2
            alpha=r2;
            if verb,     fprintf('r2 used: %g.\n',r2),end
            done=1;
        end;
    end
    
    if ~done % SW2 (and Wc's) not true at r1 or r2
        if verb,   fprintf('Must adjust bounding interval and recalc alpha. \n'),end
        if r1*r2<0 % one pos, one neg, ignore neg
            if verb,      fprintf('Only one pos point found. \n'),end
            r=max(r1,r2); % set r= the positive one
            if r==r1
                phir=phir1;
                dphir=dphir1;
            else
                phir=phir2;
                dphir=dphir2;
            end
            %Case 1
            if dphir>0
                if verb,    fprintf('Case 1; dphir >0.\n'),end
                %alpha=r; % I should check that r<=alpha before updating
                %                 phialpha=phir;
                %                 dphialpha=dphir;
                
                A=0;
                B=r;
                phiA=phi0;
                dphiA=dphi0;
                phiB=phir;
                dphiB=dphir;
                
            else %dphi<=0
                if verb,    fprintf('Case 2; dphir <=0.\n'),end
                A=r;
                B=alpha;
                phiA=phir;
                dphiA=dphir;
                phiB=phialpha;
                dphiB=dphialpha;
            end
            
            
        else % r1, r2 both positive
            if verb,   fprintf('Two  pos points found \n'),end
            if   dphir1*dphir2<0 %opposite slope condition
                if verb,     fprintf('opposite slopes.\n'),end
                A=r1;
                B=r2;
                phiA=phir1;
                dphiA=dphir1;
                phiB=phir2;
                dphiB=dphir2;
            else %dphir1*dphir2>=0 % same slope condition
                if verb,     fprintf('same slopes.\n'),end
                if dphir1<0
                    if verb,    fprintf('Both slopes neg.\n'),end
                    A=max(r1,r2);% I should check that A<=alpha before updating
                    B=alpha;
                    phiB=phialpha;
                    dphiB=dphialpha;
                    if A==r1
                        phiA=phir1;
                        dphiA=dphir1;
                    else
                        phiA=phir2;
                        dphiA=dphir2;
                    end
                    
                else
                    if verb,   fprintf('Both slopes pos.\n'),end
                    A=0;
                    B=min(r1,r2);
                    phiA=phi0;
                    dphiA=dphi0;
                    if B==r1
                        phiB=phir1;
                        dphiB=dphir1;
                    else
                        phiB=phir2;
                        dphiB=dphir2;
                    end
                end % dphir1 & dphir2 both pos
            end
        end
    end
    sw2count= sw2count+1;
end; %while not done
if sw2count>sw2max
    error('No alpha sat sw2 found.')
end

counter=counter+sw2count;

%%%%%%%%%%%%%%%% END OF TOP LEVEL CODE. FUNCTION DEFINITIONS FOLLOW %%%%%%





    function [W,W1,W2,Wapp,SW2]=WOLFEfast(alpha,c1,c2,phi0,dphi0,phialpha,dphialpha,epsk)
        W1=c1*dphi0>=(phialpha-phi0)/alpha;
        W2=dphialpha>=c2*dphi0;
        Wapp= ((2*c1-1)*dphi0>=dphialpha) && W2;
        SW2=abs(dphialpha)<=-c2*dphi0;
        W=(W1&&W2)||(Wapp &&phialpha <=phi0+epsk);
    end

    function [r1,r2,phir1,dphir1,phir2,dphir2]=invcubrat(phi,A,B,x,p,phiA,dphiA,phiB,dphiB)
        TL1=-(A*dphiA+2*A*dphiB+2*B*dphiA+B*dphiB);
        if abs(B-A)>=100*eps
            deltaphi=(phiB-phiA)/(B-A);
            BLrat = -dphiA-dphiB+2*deltaphi;
            sqrtarg=(dphiA^2+dphiA*dphiB+dphiB^2)+3*deltaphi*(-2*(dphiA+dphiB)+3*deltaphi) ;
            TL2=3*(A+B)*deltaphi;
            if abs(BLrat) >100*eps
                r1=(TL1+TL2 +  abs(B-A)*sqrt(sqrtarg ))/BLrat/3;
                r2=(TL1+TL2-   abs(B-A)*sqrt(sqrtarg ))/BLrat/3;
            else
                r1=(A+B)/2;
                r2=(A+B)/2;
            end;
        else
            r1=(A+B)/2;
            r2=(A+B)/2;
        end
        
        if nargout>2
            phir1=[];
            dphir1=[];
            phir2=[];
            dphir2=[];
            if r1>0
                [phir1,dphir1]=phi(r1,x,p);
            end
            if r2>0
                [phir2,dphir2]=phi(r2,x,p);
            end
        end
        
    end

    




    function [W,W1,W2,Wapp,SW2]=wolfe(phi, c1,c2,alpha,epsk)
        [phi0,dphi0]=phi(0,x,p);
        [phialpha,dphialpha]=phi(alpha,x,p);
        W1=c1*dphi0>=(phi(alpha,x,p)-phi0)/alpha;
        W2=dphialpha>=c2*dphi0;
        Wapp= ((2*c1-1)*dphi0>=dphialpha) && W2;
        SW2=abs(dphialpha)<=-c2*dphi0;
        W=(W1&&W2)||(Wapp &&phialpha <=phi0+epsk);
    end




    function c=init(k,x,QuadStep,alpha)
        % INIT
        % finds init vaue for c.
        
        psi0=.01;%2;
        psi1=0.1;%0.5;
        psi2=2;
        if k==0
            [f0,g0]=obj(x);
            if norm (x)>eps
                c=psi0*norm(x,'inf')/norm(g0,'inf')^2;
                return
            end
            if abs(f0)>eps
                c=psi0*abs(f0)/norm(g0,'inf')^2;
                return
            end
            c=1;
            return
        else % k>1
            [phi0,dphi0]=phi(0,x,p);
            abar=psi1*alpha;
            [dummy,dphiabar]=phi(abar,x,p);
            SConvex=(dphiabar-dphi0)>0; % quad approx to phi @ 0 & abar has a pos min
            alphaq=-dphi0*abar/(dphiabar-dphi0); % min of quad approx to phi @ 0 & abar
            if QuadStep && (phi(psi1*alpha,x,p) <=phi0) && SConvex
                c=alphaq;
                return
            else
                c=psi2*alpha;
                return
            end
        end
    end



    function [abar,bbar]=cg_update(phi,a,b,c,epsk)
        % CG_UPDATE
        % Updates (contracts) interval containing alpha.
        
        t=0.5; %Bisection method
        l=min(a,b);
        m=max(a,b);
        phi0=phi(0,x,p);
        [phic,dphic]=phi(c,x,p);
        if (c <=l) || (c>=m)
            abar= a;
            bbar = b;
            return
        end
        if  dphic >=0
            abar=a;
            bbar=c;
            return
        end
        if (phic <=phi0+epsk) %% dphic <0
            abar=c;
            bbar=b;
            return
        else %% (phic > phi0+epsk)
            abar=a;
            bbar=c;
            done=0;
            update_counter=0;
            while ~done
                done=abs(abar-bbar)<epsk;
                update_counter= update_counter+1;
                d=(1-t)*abar+t*bbar;
                [phid,dphid]=phi(d,x,p);
                if dphid >=0
                    bbar=d;
                    return
                end
                if  (phid <= phi0+epsk)
                    abar=d;
                else
                    bbar=d;
                end
            end
        end
    end



    function [abar,bbar]=cg_update_alt(phi,c,epsk)
        % CG_UPDATE_ALT
        % Alt version of Update used in CG_UPDATE_ALT.
        
        t=0.5; %Bisection method
        phi0=phi(0,x,p);
        abar=0;
        bbar=c;
        done=0;
        while ~done
            d=(1-t)*abar+t*bbar;
            [phid,dphid]=phi(d,x,p);
            if dphid >=0
                bbar=d;
                return
            end
            if  (phid <= phi0+epsk)
                abar=d;
            else
                bbar=d;
            end
        end
    end


    function [abar,bbar]=secant2(phi,a,b,epsk)
        % SECANT2
        % used to get improved estimate of alpha.
        
        
        
        c = secant (phi, a, b);
        [A, B] = cg_update (phi,a, b, c,epsk);
        if c == B,
            cbar = secant (phi,b, B);
        end
        if c == A,
            cbar = secant (phi,a, A);
        end
        if c == A || c == B,
            [abar, bbar] = cg_update (phi,A, B, cbar,epsk);
        else
            abar=A; bbar=B;
        end;
    end


    function c=secant(phi,a,b)
        % SECANT
        % Secant method to improve estimate of alpha.
        
        [f,dphia]=phi(a,x,p);
        [f,dphib]=phi(b,x,p);
        if abs(b-a)<10*eps || (abs(dphib-dphia)<10*eps)
            c=(a+b)/2;
        else
            c=(a*dphib-b*dphia)/(dphib-dphia);
        end
    end


    function [a,b]=bracket(phi,c,rho,epsk)
        % BRACKET
        % Bracket generates initial interval containing alpha.
        verb=0;
        [phi0,dphi0]=phi(0,x,p);
        j=0;
        i=-1;
        cj=c;
        ci=0; % as phici <=phi0+epsk for ci=0
        a=0;
        b=c;
        done=0;
        while ~done
            [phicj,dphicj]=phi(cj,x,p);
            if verb, fprintf('j: %i,phicj: %g, dphicj: %g \n',j,phicj,dphicj),end
            if dphicj >=0
                b=cj;
                a=ci;
                if verb,      fprintf('Exit@B1 \n'),end
                return
            end
            if phicj <=phi0+epsk
                ci=cj;
                % i=i+1;
                if verb,       fprintf('ci updated: %g, \n', ci),end
            end;
            %dphic <0
            if phicj >phi0+epsk
                [a,b]=cg_update_alt(phi,cj,epsk);
                if verb,        fprintf('Exit@B2 \n'),end
                %                pause
                return
            else %dphic <0 and  phic <= phi0+epsk
                cj=rho*cj;
                j=j+1;
                if verb,       fprintf('Exit conditions not met, inc j: %i\n',j),end
            end
        end
    end
end


