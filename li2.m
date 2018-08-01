function out=li2(y)
    x = linspace(eps,1-eps,1e6);
    dx=x(2)-x(1);
    dilog = -cumsum(log(1-x)./x)*dx;
    out=interp1(x,dilog,1-y);
end
