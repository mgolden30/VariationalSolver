function p = phase(z,params)
  [omega, ~, ~, ~] = unpack_state(z, params);

  N = params.N;
  M = params.M;

  k = 0:N-1;
  k(k>N/2) = k(k>N/2) - N;

  l = 0:M-1;
  l(l>M/2) = l(l>M/2) - M;

  o2 = fftn(omega);

  [kx, ky] = ndgrid(k);

  %reshape macros
  r1 = @(x) reshape(x, [N*N,1] );
  r2 = @(x) reshape(x, [N*N,M] );
  r3 = @(x) reshape(x, [N,N,M] );
  r4 = @(x) reshape(x, [N*N*M,1]);

  o2 = r2(o2);

  dodx = real(ifftn( r3(1i*r1(kx).*o2) ));
  dody = real(ifftn( r3(1i*r1(ky).*o2) ));
  dodt = real(ifftn( r3(1i*l     .*o2) ));

  p = [r4(dodx)'; r4(dody)'; r4(dodt)' ]; %phase matrix 
end