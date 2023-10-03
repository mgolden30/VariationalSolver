function [omega, u0, v0, T, params] = unpack_state( z, params )
  N = params.N;
  M = params.M;

  omega = reshape( z(1:M*N*N), [N,N,M] );
  u0 = z(M*N*N + 1);
  v0 = z(M*N*N + 2);
  T  = z(M*N*N + 3);
end