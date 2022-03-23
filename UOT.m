function M = UOT(eta,C,mu,nu,iter)
% Solve unbalanced OT 
% mu, nu ~measures
% C ~ loc cos fo distance between supports of mu nu
% iter ~ iterations
% eta ~ regularization parameter
K = exp(-eta*C);
[m,n] = size(K);
u = ones(m,1);
v = ones(n,1);
fi = eta/(1+eta);
for i = 1:iter
    Kv = K*v;
    u = (mu./Kv).^fi;
    Ku = K' * u;
    v = (nu./Ku).^fi;
end
M =  K .* (u * v');
end