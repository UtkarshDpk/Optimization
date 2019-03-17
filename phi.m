function [ph,dph]=phi(alpha,x,p);
%   PHI
%   phi(alpha)=f(obj(x+alpha*p)
%   phi returns phi & phi'


% if isempty(alpha)
%     error('alpha empty')
% end

if nargout==1
    f=obj(x+alpha*p);
    ph=f;
else
    [f,g]=obj(x+alpha*p);
    ph=f;
    dph=p'*g;
end
