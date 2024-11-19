
function u = controller(params, t, phi, phidot)


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

  u = kp*phi+kd*phidot +newstate*ki;
end

