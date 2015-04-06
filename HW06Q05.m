%HW06 Q05: Eigenvalue Diffusion Equation Solver

%We will define an appropriate matrix here corresponding with APhi =(1/k)FPhi
%where Phi is our Phi_i matrix and A and F have to do with the physics of the
%problem

%We will solve using the finite difference method to discretize, then 
%the power iteration to find the dominant eigenvalue,and lastly GS to finish the algorithm
%given values are a constant Sigma_a = 0.7 and
%_nu*Sigma_f = 0.6, both in cm^-1 units... and h = 0.1 cm; a = 4cm, D =1cm, 
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

hold on
plot([-a:h:a],[0;Phi;0],'ro')
%plot([-a:h:a],[0;Phi;0])
title('Eigenvalue Form of the Diffusion Equation Solution','FontSize',13)
xlabel('X,   (cm)','FontWeight','b')
ylabel('Phi/norm(Phi)','FontWeight','b')