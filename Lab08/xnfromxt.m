function xn=xnfromxt(xt,t,tn,Ts,tol)
% xt is input data
% t is times of the input data
% tn is times of the output data
% Ts is time sampling of the output data
% tol is the pinv tolerance

nt = length(t);
nn = length(tn);

A = zeros(nt,nn);

for i1=0:nt-1
    for i2=0:nn-1
        A(i1+1,i2+1) = sinc((t(i1+1)-tn(i2+1))/Ts);
    end
end

xn = pinv(A,tol)*xt;
