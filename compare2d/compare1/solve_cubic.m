function x = solve_cubic(a,b,c,d)
test = 0;
%%
x = zeros(size(b));
b = b./a;
c = c./a;
d = d./a;

p = c - b.^2/3;
q = 2*b.^3/27 - b.*c/3 + d;
delta = q.^2/4 + p.^3/27;
%%
if test == 1
    if p == 0
        x = -nthroot(q,3);
    elseif delta == 0
        x = max(3*q/p, -3*q/p/2);
    elseif delta > 0
        s = sqrt(delta);
        x = nthroot(-q/2-s,3) + nthroot(-q/2+s,3);
    else
        x = max(roots([1,0,p,q]));
%         r = 2*sqrt(-p/3);
%         s = 3*q/(p*r);
%         theta = real(acos(s)/3);
%         x = r*cos(theta);
    end
else
    %%
    ind = (p==0);
    x(ind) = -nthroot(q(ind),3);

    ind = (p~=0 & delta==0);
    x(ind) = max(3*q(ind)./p(ind), -3*q(ind)./p(ind)/2);

    ind = (p~=0 & delta>0);
    s = sqrt(delta(ind));
    x(ind) = nthroot(-q(ind)/2-s,3) + nthroot(-q(ind)/2+s,3);

    ind = (p~=0 & delta<0);
    r = 2*sqrt(-p(ind)/3);
    s = 3*q(ind)./(p(ind).*r);
    theta = real(acos(s)/3);
    x(ind) = r.*cos(theta);
end

x = x - b/3;    

end