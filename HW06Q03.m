%HW06Q03: Relative Error of the Finite Difference method for the DE

%Parameters: vary the mesh size for h = {1,.5,.1,.05,.01}
count = 1;
for h = [1,0.5,0.1,0.05,0.01]
%%%Copy HW06Q02...
a = 4; D = 1; Sigma = 0.2; S = 8; n = ((2*4)/h); L = sqrt(D/Sigma);

%lets define our A matrix, the x vector we are solving for will have Phi_1,
%Phi_2,....Phi_n-1
z = 2+(h^2/L^2);
A = zeros(n-1,n-1); b = zeros(n-1,1);

A(1,1:2) = [z, -1]; A(end, end-1:end) = [-1, z];
b(1:end,1) = (h^2)*S/D;

for i = 2:n-2
    A(i,i-1:i+1) = [-1, z, -1];
end
%now simply solve for x which we will call Phi
Phi = A\b; Phi = [0;Phi;0];

%define the function results and relative error...
q = sqrt(Sigma/D); f = @(x) (S/Sigma)*(1-((exp(q*x)+exp(-q*x))/...
    (exp(q*a)+exp(-q*a))));

Phi_actual = f([-a:h:a])';
error = abs(Phi-Phi_actual)./Phi_actual; maxerror = max(error);

compilation(count,1:2) = [h,maxerror];
count = count+1;
end
%((2*a)/h)
hold on
% set(gca,'yscale','log')
% set(gca,'xscale','log')
plot(compilation(:,1),compilation(:,2))
scatter(compilation(:,1),compilation(:,2),'ro','filled')
xlabel('Mesh Length: (h)','FontWeight','b')
ylabel('Maximum Relative Error','FontWeight','b')
title('Relative Error vs. Mesh Length',...
    'FontSize',13)
