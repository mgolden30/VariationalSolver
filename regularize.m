function z = regularize( z, params )
  %Enforce physical conditions

  N = params.N;
  M = params.M;

  [omega, u0, v0, T] = unpack_state(z, params);

  %subtract the mean
  omega = reshape( omega, [N,N,M] );
  omega = omega - mean(omega, [1,2] ); %zero mean at all times

  % enforce constant enstrophy
  %en = mean( omega(:,:,1).^2, [1,2] );
  %omega(:,:,2:M) = omega(:,:,2:M) .* sqrt(en./mean(omega(:,:,2:M).^2, [1,2]));

  %normalize so max(|omega|) = 1
  w = max(abs(omega), [], 'all' ); %vorticity scale

  %rescale via dimensional analysis
  omega = omega/w;
  u0 = u0/w;
  v0 = v0/w;
  T  = T*w;
  
  %repack
  z = pack_state( omega, u0, v0, T, params );
end