clear all;
F=@(x) (x*x)/2;
Re = input("Enter a Reynold Number either 10 or 50 : ");  % Number of grid points
disp(['Your Reynold Number : ',num2str(Re)]);
%Re = 10;           % Reynolds number (you can change this)
nu = 1 / Re;       % Viscosity
x_min = -1.0;      % Minimum x boundary
x_max = 1.0;       % Maximum x boundary
% Nx = 100;          % Number of spatial points
dx=0.08;
Nx=((x_max-x_min)/0.1)+1;
u_max=1;
dt=0.01;
t_f=1;
Nt=(t_f/dt)+1;
u=zeros(Nx,Nt);
u(1:((Nx+1)/2),1)=1;
u((((Nx+1)/2)+1):end,1)=0;
u(1,:)=1;
u(end,:)=0;
x = linspace(x_min, x_max, Nx); % Spatial grid
t = linspace(0, t_f, Nt);

u_accu=u;
u_expl=u;
u_impl=u;
u_cn=u;

% Define the exact solution function
numerator_fun1 = @(xi,x,t) ((x - xi) ./ t) .*exp(-0.5 .* Re .*(xi + (0.5 .* ((x - xi) .* (x - xi)) ./ t)));
numerator_fun2 = @(xi, x, t) ((x - xi) ./ t) .* exp(-0.5 .* Re .* 0.5 .* ((x - xi) .* (x - xi)) ./ t);
denomerator_fun1 = @(xi, x, t) exp(-0.5 .* Re .*(xi + (0.5 .* ((x - xi) .* (x - xi)) ./ t)));
denomerator_fun2 = @(xi, x, t) exp(-0.5 .* Re .* 0.5 .* ((x - xi) .* (x - xi)) ./ t);


% For Exact solution
for j=2:Nt  % time loop
    for i=2:Nx-1 %space loop
        num_integrand1=@(xi) numerator_fun1(xi,x(i),t(j));
        num_integrand2=@(xi) numerator_fun2(xi,x(i),t(j));
        den_integrand1=@(xi) denomerator_fun1(xi,x(i),t(j));
        den_integrand2=@(xi) denomerator_fun2(xi,x(i),t(j));
        numerator=integral(num_integrand1,-inf,0) + integral(num_integrand2,0,inf);
        denomerator=integral(den_integrand1,-inf,0) + integral(den_integrand2,0,inf);
        u_accu(i,j)=numerator/denomerator;
    end        
end    


% for explicit
for j=1:Nt
    for i=2:Nx-1
        u_expl(i,j+1)=u_expl(i,j) + (dt/(Re*dx*dx))*(u_expl(i-1,j) + u_expl(i+1,j) - 2*u_expl(i,j)) - (dt/dx)*(F(u_expl(i,j)) - F(u_expl(i-1,j)));
    end    
end   


%for Crank Nicolson
for j=1:Nt-1
    a=zeros(1,Nx-3);
    b=zeros(1,Nx-2);
    c=zeros(1,Nx-3);
    d=zeros(1,Nx-2);
    for i=2:Nx-2
        a(i-1)=-(0.5*dt*dx*u_cn(i,j))-(nu*dt);
        c(i-1)=(0.5*dt*dx*u_cn(i+1,j))-(nu*dt);
    end
    for i=2:Nx-1
        b(i-1)=(2*nu*dt)+(2*dx*dx);
        d(i-1)=(2*dx*dx*u_cn(i,j))+((nu*dt)*(u_cn(i+1,j)+u_cn(i-1,j)-(2*u_cn(i,j))));
    end    
    d(1)=(2*dx*dx*u_cn(2,j))+((nu*dt)*(u_cn(3,j)+u_cn(1,j)-(2*u_cn(2,j)))) + (u_cn(1,j+1)*((0.5*dt*dx*u_cn(1,j))+(nu*dt))); 
    d(Nx-2)=(2*dx*dx*u_cn(Nx-1,j))+((nu*dt)*(u_cn(Nx,j)+u_cn(Nx-2,j)-(2*u_cn(Nx-1,j)))) + (u_cn(Nx,j+1)*((nu*dt)-(0.5*dt*dx*u_cn(Nx,j))));
    x_sol=solve_tridiagonal(a,b,c,d);
    for i=2:Nx-1
        u_cn(i,j+1)=x_sol(i-1);
    end    
end    


%for Implicit
for j=1:Nt-1
    a=zeros(1,Nx-3);
    b=zeros(1,Nx-2);
    c=zeros(1,Nx-3);
    d=zeros(1,Nx-2);
    for i=2:Nx-2
        a(i-1)=-(dt*dx*u_impl(i,j))-(2*nu*dt);
        c(i-1)=(dt*dx*u_impl(i+1,j))-(2*nu*dt);
    end
    for i=2:Nx-1
        b(i-1)=(4*nu*dt)+(2*dx*dx);
        d(i-1)=(2*dx*dx*u_impl(i,j)) + (0.5*dx*dt*((u_impl(i+1,j)*u_impl(i+1,j))-(u_impl(i-1,j)*u_impl(i-1,j))));
    end  
    d(1)=(2*dx*dx*u_impl(2,j)) + (0.5*dx*dt*((u_impl(3,j)*u_impl(3,j))-(u_impl(1,j)*u_impl(1,j)))) + (u_impl(1,j+1)*((2*nu*dt)+(dt*dx)));
    d(Nx-2)=(2*dx*dx*u_impl(Nx-1,j)) + (0.5*dx*dt*((u_impl(Nx,j)*u_impl(Nx,j))-(u_impl(Nx-2,j)*u_impl(Nx-2,j)))) + (u_impl(Nx,j+1)*((2*nu*dt)-(dt*dx*u_impl(Nx,j))));
    x_sol=solve_tridiagonal(a,b,c,d);
    for i=2:Nx-1
        u_impl(i,j+1)=x_sol(i-1);
    end    
end   

figure;
hold on;
exact_line = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Exact Solution');
legend show;
xlabel('x');
ylabel('u(x,t)');
title('Moving Wave: Exact Solution');
grid on;
xlim([x_min, x_max]);
ylim([-0.1, 1.1]);

for j = 1:Nt
    clearpoints(exact_line);
    addpoints(exact_line, x, u_accu(:, j));
    % Update the plot for the current time step
    drawnow;
    % Pause to control the speed of animation (adjust if necessary)
    pause(0.05);
end

figure;
hold on;
implicit_line = animatedline('Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Implicit Solution');
legend show;
xlabel('x');
ylabel('u(x,t)');
title('Moving Wave: Implicit Solution');
grid on;
xlim([x_min, x_max]);
ylim([-0.1, 1.1]);

for j = 1:Nt
    clearpoints(implicit_line);
    addpoints(implicit_line, x, u_impl(:, j));
    % Update the plot for the current time step
    drawnow;
    % Pause to control the speed of animation (adjust if necessary)
    pause(0.05);
end

figure;
hold on;
explicit_line = animatedline('Color', 'm', 'LineWidth', 1.5, 'DisplayName', 'Explicit Solution');
legend show;
xlabel('x');
ylabel('u(x,t)');
title('Moving Wave: Explicit Solution');
grid on;
xlim([x_min, x_max]);
ylim([-0.1, 1.1]);

for j = 1:Nt
    clearpoints(explicit_line);
    addpoints(explicit_line, x, u_expl(:, j));
    % Update the plot for the current time step
    drawnow;
    % Pause to control the speed of animation (adjust if necessary)
    pause(0.05);
end

figure;
hold on;
Crank_line = animatedline('Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Crank Nicolson Solution');
legend show;
xlabel('x');
ylabel('u(x,t)');
title('Moving Wave: Crank Nicolson Solution');
grid on;
xlim([x_min, x_max]);
ylim([-0.1, 1.1]);

for j = 1:Nt
    clearpoints(Crank_line);
    addpoints(Crank_line, x, u_cn(:, j));
    % Update the plot for the current time step
    drawnow;
    % Pause to control the speed of animation (adjust if necessary)
    pause(0.05);
end


% Setup for plotting
figure;
hold on;
exact_line = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Exact Solution');
implicit_line = animatedline('Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Implicit Solution');
explicit_line = animatedline('Color', 'm', 'LineWidth', 1.5, 'DisplayName', 'Explicit Solution');
Crank_line = animatedline('Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Crank Nicolson Solution');
legend('show', 'Location', 'southwest');
xlabel('x');
ylabel('u(x,t)');
title('Moving Wave: Exact,Implicit,Explicit,Crank Nicolson');
grid on;
xlim([x_min, x_max]);
ylim([-0.1, 1.1]);

% Loop over time steps to animate the solution
for j = 1:Nt
    clearpoints(exact_line);
    clearpoints(implicit_line);
    clearpoints(explicit_line);
    clearpoints(Crank_line);
    % Add points to the exact solution line
    addpoints(exact_line, x, u_accu(:, j));
    % Add points to the implicit solution line
    addpoints(implicit_line, x, u_impl(:, j));
    addpoints(explicit_line, x, u_expl(:, j));
    addpoints(Crank_line, x, u_cn(:, j));
    % Update the plot for the current time step
    drawnow;
    % Pause to control the speed of animation (adjust if necessary)
    pause(0.01);
end

function x = solve_tridiagonal(a, b, c, d)
    n=length(d);
    c_prime=zeros(1,n-1);
    d_prime=zeros(1,n);
    c_prime(1)=c(1)/b(1);
    d_prime(1)=d(1)/b(1);
    x=zeros(1,n);

    for i = 2:n-1
        denom = b(i) - (a(i-1) * c_prime(i-1));
        c_prime(i) = c(i) / denom;
    end
    for i = 2:n
        denom = b(i) - (a(i-1) * c_prime(i-1));
        d_prime(i) = (d(i) - (a(i-1) * d_prime(i-1))) / denom;
    end

    % Backward substitution
    x(n) = d_prime(n);
    for i = n-1:-1:1
        x(i) = d_prime(i) - (c_prime(i) * x(i+1));
    end
end