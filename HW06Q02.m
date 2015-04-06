%HW06 Q02: Solving the DE using the Finite Difference Method
    %Parameters: a = 4cm, D = 1cm, Sigma = .2cm^-1, S = 8n/cm^3*s, h = .1cm
    %Equation: Phi'' - 1/L^2Phi = -S/D, L = sqrt(D/Sigma)
    %B.C.'s: Phi(a) = Phi(-a) = 0; solve on x E [-a,a]
    
a = 4; D = 1; Sigma = 0.2; S = 8; h = .1; n = ((2*a)/h); L = sqrt(D/Sigma);

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

Phi = A\b;

Phi = [0;Phi;0]; x = [-a:h:a]';

hold on
plot(x,Phi,'bo')
xlabel('X   [cm]', 'FontWeight','b'); ylabel('Phi(x)   [1/(cm^2*s)]',...
    'FontWeight','b')
title('Diffusion Equation Solution: Finite Difference Method',...
    'FontWeight','b', 'FontSize',13)

%now lets compare to HW06Q01_b b/c essentially the exact same situations...

q = sqrt(Sigma/D); f = @(x) (S/Sigma)*(1-((exp(q*x)+exp(-q*x))/...
    (exp(q*a)+exp(-q*a))));

x_analytic = [-a:h:a]; F = f(x_analytic);

plot(x_analytic,F,'r-')
legend('Finite Difference','Analytic Solution')
