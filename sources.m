function u = sources(N,omega,tol)
% sources(N,omega,tol)
% Use SOR to solve diffusion with sources and sinks
% N = number of samples
% omega = relaxation parameter 0<omega<2
% tol = max residual tolerance

global M
M = N;

lambda=100.;
x = linspace(-2,2,2*N-1);
y = linspace(-1,1,N);
[X,Y] = meshgrid(x,y);
X = X';
Y = Y';
dx = x(2)-x(1);
dy = y(2)-y(1);
g = 10*lambda/sqrt(pi)*exp(-lambda^2*((X-1).^2+Y.^2))-10*lambda/sqrt(pi)*exp(-lambda^2*((X+1).^2+Y.^2));

u = zeros(size(X));

figure(1)
subplot(2,1,1);
contourf(X,Y,u);
axis equal
subplot(2,1,2);
surf(X,Y,u);
axis equal
drawnow

maxresid = 1+tol
iter = 0;
while maxresid > tol
    maxresid = 0.;
    i = 1;
    for j=2:N-1
        resid = 0.25*(u(i+1,j)+u(i,j+1)+u(i,j-1)-3*u(i,j))-dx^2/4*g(i,j);
        maxresid = max([abs(resid),maxresid]);
        u(i,j) = u(i,j) + omega*resid;
    end
    
    for i=2:2*N-2
        for j=2:N-1
            resid = 0.25*(u(i+1,j)+u(i,j+1)+u(i,j-1)+u(i-1,j)-4*u(i,j))-dx^2/4*g(i,j);
            maxresid = max([abs(resid),maxresid]);
            u(i,j) = u(i,j) + omega*resid;
        end
    end
    
    i = 2*N-1;
    for j=2:N-1
        resid = 0.25*(u(i-1,j)+u(i,j+1)+u(i,j-1)-3*u(i,j))-dx^2/4*g(i,j);
        maxresid = max([abs(resid),maxresid]);
        u(i,j) = u(i,j) + omega*resid;
    end
    
    subplot(2,1,1);
    contourf(X,Y,u);
    axis equal
    title(sprintf('maxresid = %e',maxresid));
    subplot(2,1,2);
    surf(X,Y,u);
    drawnow
    iter = iter+1;
end

iter

end

function r = col(i,j)

global M

r = (2*M+1)*(j-1)+i;
end