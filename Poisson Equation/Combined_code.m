S = @(x, y) (2*sinh(10*(x-0.5))) + (40*(x-0.5)*cosh(10*(x-0.5))) + (100*(x-0.5).^2*sinh(10*(x-0.5))) + ...
             (2*sinh(10*(y-0.5))) + (40*(y-0.5)*cosh(10*(y-0.5))) + (100*(y-0.5).^2*sinh(10*(y-0.5))) + ...
             (4*(x.^2 + y.^2)*exp(2*x*y));

Anal_sol = @(x, y) ((x-0.5).^2*sinh(10*(x-0.5))) + ((y-0.5).^2*sinh(10*(y-0.5))) + (exp(2*x*y));

N = input("Enter a Grid Size either 41 or 81 : ");  % Number of grid points
disp(['Your Grid Size : ',num2str(N)]);
L = 1;   % Domain length
h = L / (N - 1);  % Grid spacing

% Define phi with boundary conditions
phi = ones(N, N);

% Boundary condition functions
f1 = @(y) ((y-0.5).^2*sinh(10*(y-0.5))) + (0.25*sinh(-5)) + 1;
f2 = @(y) ((y-0.5).^2*sinh(10*(y-0.5))) + (0.25*sinh(5)) + exp(2*y);
f3 = @(x) ((x-0.5).^2*sinh(10*(x-0.5))) + (0.25*sinh(-5)) + 1;
f4 = @(x) ((x-0.5).^2*sinh(10*(x-0.5))) + (0.25*sinh(5)) + exp(2*x);

% Apply boundary conditions
for i = 1:N
    phi(1, i) = f1((i-1)*h);    % Bottom boundary
    phi(N, i) = f2((i-1)*h);    % Top    boundary
    phi(i, 1) = f3((i-1)*h);    % Left   boundary
    phi(i, N) = f4((i-1)*h);    % Right  boundary
end

phi_cons=phi;
max_iter = 50000;
tol = 1e-6;

% Right hand side-matrix
f = zeros(N, N);
for j = 1:N
    for i = 1:N
        f(i, j) = S((i-1)*h, (j-1)*h);
    end
end


% Jacobi method
phi=phi_cons;
phi_new_jacobi=phi;
%x_jacobi=zeros(1,200);
l=1;
for k=1:max_iter
    phi_old=phi;
    for j=2:N-1
        for i=2:N-1
            phi_new_jacobi(i,j)=0.25*(phi_old(i+1,j)+phi_old(i-1,j)+phi_old(i,j+1)+phi_old(i,j-1)-(h^2 * f(i,j)));
        end    
    end
    if sum(sum(abs(phi_new_jacobi - phi_old)))<0.1
        residual_jacobi(l)=sum(sum(abs(phi_new_jacobi - phi_old)));
        x_jacobi(l)=k;
        l=l+1;
    end    
    if sum(sum(abs(phi_new_jacobi - phi_old))) < tol
        fprintf("Jacobi Iteration : %d\n",k);
        break;
    end
    phi=phi_new_jacobi;
end 

% Gauss Seidal
phi=phi_cons;
phi_new_gauss=phi;
%residual_gauss=zeros(1,residual_iteration);
%x_Gauss=zeros(1,200);
l=1;
for k=1:max_iter
    phi_old=phi;
    for j=2:N-1
        for i=2:N-1
            phi_new_gauss(i,j)=0.25*(phi_old(i+1,j)+phi_new_gauss(i-1,j)+phi_old(i,j+1)+phi_new_gauss(i,j-1)-(h^2 * f(i,j)));
        end    
    end
    if sum(sum(abs(phi_new_gauss - phi_old)))<0.1
        residual_gauss(l)=sum(sum(abs(phi_new_gauss - phi_old)));
        x_Gauss(l)=k;
        l=l+1;
    end  
    if sum(sum(abs(phi_new_gauss - phi))) < tol
        fprintf("Gauss Seidal Iteration : %d\n",k);
        break;
    end
    phi=phi_new_gauss;
end 



% Line SOR parameters
a = 1;
b = -4;
c = 1;
rel_fac = 1.2;  % Start with relaxation factor 1.0 and adjust later if necessary
phi=phi_cons;
phi_new_LSOR = phi;
phi_new_LSOR1 = phi;
%residual_LSOR=zeros(1,residual_iteration);
%x_LSOR=zeros(1,200);
l=1;

% Start Line SOR iterations
for k = 1:max_iter
    phi_old = phi;  % Store old solution for comparison later
    
    % Solve row by row using tridiagonal solver
    for j = 2:N-1
        d = zeros(1, N-2);
        for i = 2:N-1
            d(i-1) = (h*h*f(i,j)) - (phi_new_LSOR(i,j-1)) - phi_old(i,j+1);
        end
        d(1)=(h*h*f(2,j))-(phi_new_LSOR(2,j-1))-(phi_old(2,j+1))-(phi_new_LSOR(1,j));
        d(N-2)=(h*h*f(N-1,j))-(phi_new_LSOR(N-1,j-1))-(phi_old(N-1,j+1))-(phi_new_LSOR(N,j));

        % Solve tridiagonal system
        x = solve_tridiagonal(a, b, c, d, N-2);
        
        % Update phi_new with the solution from tridiagonal solver
        for i = 2:N-1
            phi_new_LSOR(i,j) = x(i-1);
        end
    end
    
    % Apply relaxation factor
    for j = 2:N-1
        for i = 2:N-1
            phi_new_LSOR1(i,j) = rel_fac * phi_new_LSOR(i,j) + (1-rel_fac) * phi_old(i,j);
        end
    end
    
    %Check convergence
    
    if sum(sum(abs(phi_new_LSOR1 - phi_old)))<0.1
        residual_LSOR(l)=sum(sum(abs(phi_new_LSOR1 - phi_old)));
        x_LSOR(l)=k;
        l=l+1;
    end     
    if sum(sum(abs(phi_new_LSOR1 - phi_old))) < tol
        fprintf("LSOR Iteration : %d\n",k);
        break;
    end
    
    % Update phi for the next iteration
    phi = phi_new_LSOR1;
end

% RBPGS parameters
phi_RBPGS=phi_cons;
phi_old_RBPGS=phi_cons;
%residual_RBPGS=zeros(1,residual_iteration);
%x_RBPGS=zeros(1,200);
l=1;

for k=1:max_iter
    for j=2:N-1
        for i=2:N-1
            if(rem((i+j),2)==0)
                phi_RBPGS(i,j)=(phi_old_RBPGS(i-1,j)+phi_old_RBPGS(i+1,j)+phi_old_RBPGS(i,j-1)+phi_old_RBPGS(i,j+1)-(h^2*f(i,j)))/4;
            end    
        end    
    end

    for j=2:N-1
        for i=2:N-1
            if(rem((i+j),2)==1)
                phi_RBPGS(i,j)=(phi_RBPGS(i-1,j)+phi_RBPGS(i+1,j)+phi_RBPGS(i,j-1)+phi_RBPGS(i,j+1)-(h^2*f(i,j)))/4;
            end    
        end    
    end

    if sum(sum(abs(phi_RBPGS - phi_old_RBPGS)))<0.1
        residual_RBPGS(l)=sum(sum(abs(phi_RBPGS - phi_old_RBPGS)));
        x_RBPGS(l)=k;
        l=l+1;
    end  

    if sum(sum(abs(phi_RBPGS - phi_old_RBPGS))) < tol
        fprintf("RBPGS Iteration : %d\n",k);
        break;
    end
    phi_old_RBPGS=phi_RBPGS;
end


% ADI parameters
a = 1;
b = -4;
c = 1;
phi=phi_cons;
phi_new_ADI = phi;
phi_new_ADI_1 = phi;
%residual_ADI=zeros(1,residual_iteration);
%x_ADI=zeros(1,200);
l=1;

% Start ADI iterations
for k = 1:max_iter
    phi_old = phi;  % Store old solution for comparison later
    
    % First half solve in row wise x direction 
    for j = 2:N-1
        d = zeros(1, N-2);
        for i = 2:N-1
            d(i-1) = (h*h*f(i,j)) - (phi_new_ADI(i,j-1)) - phi_old(i,j+1);
        end
        d(1)=(h*h*f(2,j))-(phi_new_ADI(2,j-1))-(phi_old(2,j+1))-(phi_new_ADI(1,j));
        d(N-2)=(h*h*f(N-1,j))-(phi_new_ADI(N-1,j-1))-(phi_old(N-1,j+1))-(phi_new_ADI(N,j));

        % Solve tridiagonal system
        phi_new_row = solve_tridiagonal(a, b, c, d, N-2);
        
        % Update phi_new with the solution from tridiagonal solver
        for i = 2:N-1
            phi_new_ADI(i,j) = phi_new_row(i-1);
        end
    end
    
    % For second half column sweep
    for i=2:N-1
        d=zeros(1,N-2);
        for j=2:N-2
            d(j-1)=(h*h*f(i,j))-(phi_new_ADI(i+1,j))-(phi_new_ADI_1(i-1,j));
        end
        d(1)=(h*h*f(i,2))-(phi_new_ADI_1(i-1,2))-(phi_new_ADI(i+1,2))-(phi_new_ADI_1(i,1));
        d(N-2)=(h*h*f(i,N-1))-(phi_new_ADI_1(i-1,N-1))-(phi_new_ADI(i+1,N-1))-(phi_new_ADI(i,N));

        %Solve for tridiagonal system
        phi_new_col=solve_tridiagonal(a,b,c,d,N-2);

        % Update phi_new with the solution from tridiagonal solver
        for j = 2:N-1
            phi_new_ADI_1(i,j)=phi_new_col(j-1);
        end    

    end    

    
    %Check convergence
    if sum(sum(abs(phi_new_ADI_1 - phi_old)))<0.1
        residual_ADI(l)=sum(sum(abs(phi_new_ADI_1 - phi_old)));
        x_ADI(l)=k;
        l=l+1;
    end     
    if sum(sum(abs(phi_new_ADI_1 - phi_old))) < tol
        fprintf("ADI Iteration : %d\n",k);
        break;
    end
    
    % Update phi for the next iteration
    phi = phi_new_ADI_1;
    
end



[X, Y] = meshgrid(0:h:L, 0:h:L);

% Exact solution computed using the analytical function S
exact_solution = zeros(N, N);
for j = 1:N
    for i = 1:N
        exact_solution(i, j) = Anal_sol((i-1)*h, (j-1)*h);
    end
end

% Plot both exact and numerical solutions on the same figure
%Figure 1;
% Figure 1: Surface plots comparing exact and numerical solutions
figure(1);

% Plot exact solution
subplot(3, 2, 1);
surf(X, Y, exact_solution);
title('Exact Solution');
xlabel('x');
ylabel('y');
zlabel('phi(x, y)');
colorbar;

% Plot numerical solution of Jacobi
subplot(3, 2, 2);
surf(X, Y, phi_new_jacobi);
title('Numerical Solution (Jacobi)');
xlabel('x');
ylabel('y');
zlabel('phi(x, y)');
colorbar;

% Plot numerical solution of Gauss-Seidel
subplot(3, 2, 3);
surf(X, Y, phi_new_gauss);
title('Numerical Solution (Gauss-Seidel)');
xlabel('x');
ylabel('y');
zlabel('phi(x, y)');
colorbar;

% Plot numerical solution of Line SOR
subplot(3, 2, 4);
surf(X, Y, phi_new_LSOR1);
title('Numerical Solution (Line SOR)');
xlabel('x');
ylabel('y');
zlabel('phi(x, y)');
colorbar;

% Plot numerical solution of RBPGS
subplot(3, 2, 5);
surf(X, Y, phi_RBPGS);
title('Numerical Solution (RBPGS)');
xlabel('x');
ylabel('y');
zlabel('phi(x, y)');
colorbar;

% Plot numerical solution of ADI
subplot(3, 2, 6);
surf(X, Y, phi_new_ADI_1);
title('Numerical Solution (ADI)');
xlabel('x');
ylabel('y');
zlabel('phi(x, y)');
colorbar;

% Figure 2: Line plot comparing residuals
figure(2);
%k = linspace(1, residual_iteration, residual_iteration);
plot(x_jacobi, residual_jacobi, '-r', 'LineWidth', 2);  
hold on;                           
plot(x_Gauss, residual_gauss, '--b', 'LineWidth', 2);  
plot(x_LSOR, residual_LSOR, ':g', 'LineWidth', 2);  
plot(x_RBPGS, residual_RBPGS, '-m', 'LineWidth', 2); 
plot(x_ADI, residual_ADI, '-.k', 'LineWidth', 2); 
hold off;                 

% Add labels, title, and legend
xlabel('Iterations');
ylabel('Residual');
title('Comparison of Residuals');
legend('Jacobi', 'Gauss-Seidel', 'Line SOR', 'RBPGS', 'ADI');
grid on;

% Figure 3: Contour plots comparing exact and numerical solutions
figure(3);

% Plot contour for Jacobi comparison
subplot(3, 2, 1);
contourf(X, Y, exact_solution); 
% hold on;
% contour(X, Y, phi_new_jacobi);
title('Exact Solution');
%legend('Exact Solution', 'Numerical Solution');
xlabel('x');
ylabel('y');
colorbar;

% Plot contour for Gauss-Seidel comparison
subplot(3, 2, 2);
% contour(X, Y, exact_solution, 'LineWidth', 2); 
% hold on;
contourf(X, Y, phi_new_jacobi);
title('Jacobi Solution');
%legend('Exact Solution', 'Numerical Solution');
xlabel('x');
ylabel('y');
colorbar;

% Plot contour for LSOR comparison
subplot(3, 2, 3);
% contour(X, Y, exact_solution, 'LineWidth', 2); 
% hold on;
contourf(X, Y, phi_new_gauss);
title('Gauss Seidal Solution');
%legend('Exact Solution', 'Numerical Solution');
xlabel('x');
ylabel('y');
colorbar;

% Plot contour for RBPGS comparison
subplot(3, 2, 4);
% contour(X, Y, exact_solution, 'LineWidth', 2); 
% hold on;
contourf(X, Y, phi_new_LSOR1);
title('LSOR Solution');
%legend('Exact Solution', 'Numerical Solution');
xlabel('x');
ylabel('y');
colorbar;

% Plot contour for ADI comparison
subplot(3, 2, 5);
% contour(X, Y, exact_solution, 'LineWidth', 2); 
% hold on;
contourf(X, Y, phi_RBPGS);
title('RBPGS Solution');
% legend('Exact Solution', 'Numerical Solution');
xlabel('x');
ylabel('y');
colorbar;


% Plot contour of numerical solution (RBPGS)
subplot(3, 2, 6);
% contour(X, Y, exact_solution, 'LineWidth', 2); 
% hold on;
contourf(X, Y, phi_new_ADI_1);
title('ADI Solution');
% legend('Exact Solution', 'Numerical Solution');
xlabel('x');
ylabel('y');
colorbar;



% Function to solve a tridiagonal system using the Thomas algorithm
function x = solve_tridiagonal(a, b, c, d, n)
    c_prime = zeros(1, n-1);
    d_prime = zeros(1, n);
    x = zeros(1, n);
    
    % Forward sweep
    c_prime(1) = c / b;
    d_prime(1) = d(1) / b;
    
    for i = 2:n-1
        denom = b - (a * c_prime(i-1));
        c_prime(i) = c / denom;
    end
    for i = 2:n
        denom = b - (a * c_prime(i-1));
        d_prime(i) = (d(i) - (a * d_prime(i-1))) / denom;
    end

    % Backward substitution
    x(n) = d_prime(n);
    for i = n-1:-1:1
        x(i) = d_prime(i) - (c_prime(i) * x(i+1));
    end
end
