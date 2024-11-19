
function u = controllerNoisy(params, t, obs)

  % Now you only receive noisy measurements for phi
  % obs = [ay; az; gx] with a* in units of g's, and gx in units of rad/s

  % This template code calls the function EKFupdate that you must complete below
  xhat = EKFupdate(params, t, obs);
  phi = xhat(1);
  phidot = xhat(2);

   
  persistent newstate t_last
  if isempty(newstate)
     t_last=0;
     newstate = 0;
  end
  dt=t-t_last;
  newstate=newstate+phi*dt;
  t_last=t;

  ki=9000;
  kd=0.3;
  kp =90;
 % kd=5;
  u = kp*phi+kd*phidot +newstate*ki;
end

function xhatOut = EKFupdate(params, t, z)
  % z = [ay; az; gx] with a* in units of g's, and gx in units of rad/s
 
   
        persistent t_last xhat P
        if isempty(P)
            t_last=0;
            xhat=[0;0];
            P=eye(2);
        end
        dt=t-t_last;
        
 
        A=[1 dt;0 1];
        xhat=A*xhat;
        
        P=A*P*A.';
        phi=xhat(1);
        phi_dot=xhat(2);
        H = [cos(phi) -sin(phi) 0; 0 0 1]'; 
        h = [sin(phi);cos(phi);phi_dot];
        
        K = P*H.'/(H*P*H.');
        xhat=  xhat+ K * (z - h);
        P = (eye(2) - K*H) * P; 
        
        t_last=t;
        xhatOut=xhat;
end
