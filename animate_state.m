function animate_state( z, params )
   for t = 1:params.M
      %tiledlayout(1,2);
      
      %nexttile
      [omega, ~, ~, ~, ~] = unpack_state( z, params );
      imagesc( omega(:,:,t) );
      clim([-1 1]);
      colorbar;
      set(gca, 'ydir', 'normal' );
      colormap bluewhitered;
      axis square
      title( "" + t );

      %{
      nexttile
      o = reshape(omega, [N*N,M]);

      dist = vecnorm(o - o(:,1));
      dist = [dist, dist(1)];
      plot(dist, 'o');
      %}
      
      drawnow
    end
end