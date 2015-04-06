%HW06Q06

%from question 5.....
nuxSigma_f = 0.6; Sigma_a = 0.7; h = 0.1; a = 4; D = 1; n = (2*a)/h;

z = (2*D)/(h^2) + Sigma_a; q = D/(h^2);
%Finite Difference, must define A, and F
F = zeros(n-1,n-1)+(nuxSigma_f*eye(n-1));
A = zeros(n-1,n-1);
A(1,1:2) = [z,-q]; A(end,end-1:end) = [-q,z];
for i = 2:n-2
    A(i,i-1:i+1) = [-q,z,-q];
end

%now for our Power iteration method and GS...
k_before = 1; Phi_before = ones(n-1,1);
Phi_before = Phi_before/norm(Phi_before); %normalize initial guess
Q_fbefore = F*Phi_before;

%do our first iteration to be able to get an error and enter the loop...
%Gauss-Seidel Iteration Solver: A = L + U + D; (D+L)x_k = -Ux_k-1 + b;
%P = -(D+L)^-1*U
b = (1/k_before)*Q_fbefore;
D = diag(diag(A)); D_inv = inv(D); L = tril(A) - D; U = triu(A) - D;
P = -inv(D+L)*U; b_tilda = inv(D+L)*b;

Phi_next = (P*Phi_before)+(b_tilda);
Q_fnext = F*Phi_next;
k_next = k_before*(sum(Q_fnext)/sum(Q_fbefore)); b_tilda = inv(D+L)*(1/k_next)*Q_fnext;
abserrorphi = norm(abs(Phi_next - Phi_before)); abserrork = abs(k_next-k_before);
count = 1;
k_before = k_next; Q_fbefore = Q_fnext;
%now do the main iteration part
while abserrork > 1e-4 | abserrorphi > 1e-4
    Phi_final = (P*Phi_next)+(b_tilda);
    abserrorphi = norm(abs(Phi_final - Phi_next));
    Phi_next = Phi_final; Q_fnext = F*Phi_next; 
    k_next = k_before*(sum(Q_fnext)/sum(Q_fbefore));
    b_tilda = inv(D+L)*((1/k_next)*Q_fnext);
    abserrork = abs(k_next-k_before);
    count = count+1;
    
    k_before = k_next; Q_fbefore = Q_fnext;
end
iterations = count; k = k_before; Phi = Phi_next/norm(Phi_next);
B = inv(A)*F;
%now to continue on with this question 6....
PhiNew = inv(B)*max(eig(B))*Phi;

norm_maxrelative = norm(abs((PhiNew-Phi)./(PhiNew)));
norm_absolute = norm(abs(PhiNew-Phi));

hold on
plot([-a:h:a],[0;Phi;0],'k','LineWidth',2); plot([-a:h:a],[0;PhiNew;0],'g','LineWidth',2)
xlabel('X,   (cm)','FontWeight','b')
ylabel('Phi/norm(Phi)','FontWeight','b')
title('Matlab Eigenvalue Solution vs. Numerical E-Value Solution','FontSize',13)
legend('Phi Numerical','Matlab Phi')