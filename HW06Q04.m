%HW06Q04: Iterative Solutions of the DE

%Note: we will be using the matrix created in HW06Q02, we will generate the
%number of iterations necessary for the Jacobi, Gauss-Seidel, and Succesive
%Over Relaxation for two errors, {1E-3, 1E-5} and four different mesh
%sizes, h = {1,0.5,0.1,0.05}
index = 1;
for error = [1e-3,1e-5]
    
for h = [1,0.5,0.1,0.05]
    %%%Copy HW06Q02...
a = 4; D = 1; Sigma = 0.2; S = 8; n = ((2*4)/h); L = sqrt(D/Sigma);
x_0 = zeros(n-1,1);
%lets define our A matrix, the x vector we are solving for will have Phi_1,
%Phi_2,....Phi_n-1
z = 2+(h^2/L^2);
A = zeros(n-1,n-1); b = zeros(n-1,1);

A(1,1:2) = [z, -1]; A(end, end-1:end) = [-1, z];
b(1:end,1) = (h^2)*S/D;

for i = 2:n-2
    A(i,i-1:i+1) = [-1, z, -1];
end

%Jacobi Iteration Solver: P = D^-1(D - A) = I - D^-1*A; D = diag(A)
    %we will simply define our needed matrices, with the appropriate
    %relationships provided above and then execute the iteration methods
D = diag(diag(A)); D_inv = inv(D); P = D_inv*(D-A); b_tilda = D_inv*b;

x_k = x_0;
x_knext = (P*x_0)+(b_tilda);
abserror = norm(abs(x_knext - x_k));
count = 0;
while abserror > error
    x_knext = (P*x_k)+(b_tilda);
    abserror = norm(abs(x_knext - x_k));
    x_k = x_knext;
    count = count+1;
end
x_J = x_knext; iterations_J = count;

%Gauss-Seidel Iteration Solver: A = L + U + D; (D+L)x_k = -Ux_k-1 + b;
%P = -(D+L)^-1*U
L = tril(A) - D; U = triu(A) - D;
P = -inv(D+L)*U; b_tilda = inv(D+L)*b;

x_k = x_0;
x_knext = (P*x_0)+(b_tilda);
abserror = norm(abs(x_knext - x_k));
count = 0;
while abserror > error
    x_knext = (P*x_k)+(b_tilda);
    abserror = norm(abs(x_knext - x_k));
    x_k = x_knext;
    count = count+1;
end

x_GS = x_knext; iterations_GS = count;

%Successive Over Relaxation Iteration Solver: let's multiply through by
%omega (0<w<2) to try and speed up our method a bit...
%(D+wL)x = ((1-w)D-wU)x + wb  -> P = (D+wL)^-1*((1-w)D-wU)
w = 1.2; P = inv(D+w*L)*((1-w)*D-(w*U)); b_tilda = inv(D+w*L)*w*b;

x_k = x_0;
x_knext = (P*x_0)+(b_tilda);
abserror = norm(abs(x_knext - x_k));
count = 0;
while abserror > error
    x_knext = (P*x_k)+(b_tilda);
    abserror = norm(abs(x_knext - x_k));
    x_k = x_knext;
    count = count+1;
end

x_SOR = x_knext; iterations_SOR = count;

MATRIX(index,1:5) = [error,h,iterations_J,iterations_GS,iterations_SOR];
index = index+1;
end

end
%Note, MATRIX is ordered like so: in column 1 is the step size, in column 2
%the error tolerance is listed, in column 3 the Jacobi iteration numbers,
%in column 4 are the GS iteration counts, and in Column 5 the SOR counts;
%the 8 rows are the 4 mesh lengths repeated twice for the two different
%error tolerances
MATRIX, LW = 'LineWidth';
hold on
subplot(2,1,1)
title('Iterative Solver Speed vs. Mesh Length [Error = 10e-3]','FontSize',13)
ylabel('Iterations','FontWeight','b'); xlabel('Mesh Length','FontWeight','b')
set(gca,'yscale','log')
hold on
plot([1,0.5,0.1,0.05],MATRIX(1:4,3),'r',LW,2);plot([1,0.5,0.1,0.05],MATRIX(1:4,4),'b',LW,2);
plot([1,0.5,0.1,0.05],MATRIX(1:4,5),'c',LW,2);
legend('Jacobi','Gauss-Seidel','SOR')
scatter([1,0.5,0.1,0.05],MATRIX(1:4,3),'ro','filled');scatter([1,0.5,0.1,0.05],MATRIX(1:4,4),'bo','filled');
scatter([1,0.5,0.1,0.05],MATRIX(1:4,5),'co','filled');
hold off
subplot(2,1,2)
title('Iterative Solver Speed vs. Mesh Length [Error = 10e-5]','FontSize',13)
ylabel('Iterations','FontWeight','b'); xlabel('Mesh Length','FontWeight','b')
set(gca,'yscale','log')
hold on
plot([1,0.5,0.1,0.05],MATRIX(5:8,3),'r',LW,2);plot([1,0.5,0.1,0.05],MATRIX(5:8,4),'b',LW,2);
plot([1,0.5,0.1,0.05],MATRIX(5:8,5),'c',LW,2);
legend('Jacobi','Gauss-Seidel','SOR')
scatter([1,0.5,0.1,0.05],MATRIX(5:8,3),'ro','filled');scatter([1,0.5,0.1,0.05],MATRIX(5:8,4),'bo','filled');
scatter([1,0.5,0.1,0.05],MATRIX(5:8,5),'co','filled');
hold off