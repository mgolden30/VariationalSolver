function f = ECS_objective( z, params )
  [omega, u0, v0, T] = unpack_state(z, params);

  N = params.N;
  M = params.M;

  k = 0:N-1;
  k(k>N/2) = k(k>N/2) - N;

  l = 0:M-1;
  l(l>M/2) = l(l>M/2) - M;

  o2 = fft2(omega);

  [kx, ky] = ndgrid(k);

  %define operators to compute velocities
  to_u = -1i*ky./(k.^2 + k.'.^2); % \partial_y psi
  to_v =  1i*kx./(k.^2 + k.'.^2); %-\partial_x psi

  to_u(1,1) = 0; %no mean flow
  to_v(1,1) = 0; %no mean flow

  %reshape macros
  r1 = @(x) reshape(x, [N*N,1] );
  r2 = @(x) reshape(x, [N*N,M] );
  r3 = @(x) reshape(x, [N,N,M] );
  r4 = @(x) reshape(x, [N*N*M,1]);

  u = r1(to_u).* r2( o2 );
  v = r1(to_v).* r2( o2 );

  u = reshape( u, [N,N,M] );
  v = reshape( v, [N,N,M] );

  u = real(ifft2(u)) + u0;
  v = real(ifft2(v)) + v0;

  %Now compute advection
  uw2 = fft2( u.*omega );
  vw2 = fft2( v.*omega );

  rhs = 1i*(kx.*uw2 + ky.*vw2);
  %integrate rhs in time

  time_int = 1./(1i*l);
  time_int(1) = 1; %preserve the mean

  rhs = r2(rhs);
  rhs = ifft( time_int.*fft(rhs,[],2), [], 2 );

  f2 = 2*pi/T*( r2(o2) - mean(r2(o2), 2) ) + rhs; 
  
  A = 1./( k.^2 + k.'.^2);
  A(1,1) = 1;

  A = r1(A);

  f2 = A.*f2; %take inverse Laplacian

  f2 = r3(f2); %reshape it

  f = real(ifft2(f2)); %turn it back to spatial

  f = r4(f); %reshape it
end