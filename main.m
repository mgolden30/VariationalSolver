%{
Variational solver for Periodic Orbits of the Euler equations
%}

clear;

N = 128; %points in space
M = 32; %points in time

%1d grids
space = (0:(N-1))*2*pi/N;
time  = (0:(M-1))*2*pi/M;

%3d coordinate matrices
[x,y,t] = ndgrid( space, space, time );

%Guess of periodic orbit
omega = cos(x).*cos(t) + cos(x+y).*sin(t+1);

% allowing mean flow is equivalent to continuous drift in space
% Therefore without loss of generality, all ECS are periodic orbits rather
% than RPOs
u0 =  0.1;
v0 = -0.2;

%params is just used to track the size of the fields in space and time
params.N = N;
params.M = M;

%guess at period
T = 10;

z = pack_state( omega, u0, v0, T, params );

%Newton parameters
maxit = 64;
hook  = 0.5;

%GMRES parameters
inner = 64;
outer = 1;
tol   = 1e-6;
every = 16; %how often to animate state

err = zeros(maxit,1);

for i = 1:maxit

   %check if we need to animate the trajectory
  if( mod(i-1,every) == 0 )
    %enforce constraints
    z = regularize( z, params );
    %animate_state(  z, params )
  end

  % compute two versions of objective
  % FF has no phase conditions, F does
  FF= @(z) [ECS_objective( z, params );];
  F = @(z) [FF(z); 0; 0; 0]; %three phase constraints

  %actually evaluate the functions
  ff= FF(z);
  f = F(z);
  
  %compute magnitude to monitor Newton
  err(i) = norm(ff);

  h = 1e-6; %finite difference parameter
  p = phase(z,params); %compute phase only once

  %finite differenced Jacobian + phase conditions
  J = @(v) [(FF(z+h*v) - ff)/h; p * v(1:end-3)];

  %initial guess of GMRES solution
  dz0 = 0*z;
  
  %Do GMRES
  [dz, ~] = gmres( J, f, inner, tol, outer, [], [], dz0 );

  %Print status
  fprintf("%d: |F| = %.3e, GMRES residual = %.3e\n", i, err(i)/norm(z(1:M*N*N)), norm(J(dz) - f)/err(i) );

  %Take Newton step!
  z = z - hook*dz;
end


%%
clf

tiledlayout(1,2);
nexttile
semilogy( 1:maxit, err, 'o' );

nexttile
animate_state(  z, params );
