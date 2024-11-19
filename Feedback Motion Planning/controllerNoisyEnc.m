
function u = controllerNoisyEnc(params, t, obs, th, dth)
  % Now you only receive noisy measurements for theta, 
  % obs = [ay; az; gx] 
  % New for 6b: you also have access to params.traj(t)


  xhat = EKFupdate(params, t, obs);
  phi = xhat(1);
  dphi = xhat(2);
  x=params.r*(phi+th);
  xd=params.r*(dphi+dth);
  
  kp_x=0.2;
  kd_x=0.5;

  e_x=params.traj(t)-x;
  e_xdot=0-xd;
  u_x=kp_x*(e_x)+kd_x*(e_xdot);
  phides=asin(u_x);
  
  kp_phi=0.1;
  kd_phi=0.01;
  e_phi=phides- phi ;
  e_phidot=0-dphi;
  
  uphi=kp_phi*sin(e_phi)+kd_phi*(e_phidot);
  u=-uphi;





end

function xhatOut = EKFupdate(params, t, z)
  % z = [ay; az; gx] with a* in units of g's, and gx in units of rad/s
  % using persistent variables to create/update any additional state that you need.

  Q = diag([0.01 0.01]);
  R = diag([0.1 0.1 0.1]);
  
  persistent xhat P t_ekf
  if isempty(P)
      P = eye(2);
      t_ekf = 0;
      xhat = [0;0];
  end
  
  dt = t - t_ekf;
  A = [1 dt; 0 1];
  xhat = A*xhat;
  P = A*P*A' + Q;
  
  phi = xhat(1);
  phidot = xhat(2);
  
  H = [cos(phi) -sin(phi) 0; 0 0 1]'; 
  K = P*H'/(H*P*H' + R);
  
  h = [sin(phi) cos(phi) phidot]';
  xhat = xhat + K*(z - h);
  P = (eye(2) - K*H)*P;
  
  t_ekf = t;
  xhatOut = xhat;
end
