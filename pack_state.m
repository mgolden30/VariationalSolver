function z = pack_state( omega, u0, v0, T, params)
  M = params.M;
  N = params.N;
  z = [ reshape(omega, [N*N*M,1]); u0; v0; T ];
end