
function u = controller(params, t, X)
  % You have full state feedback available
  params = struct();

  params.g = 9.81;
  params.mr = 0.25;
  params.ir = 0.0001;
  params.d = 0.1;
  params.r = 0.02;

  % you will need the state variables, and control input to be declared as symbolic
  syms th phi dth dphi u


  qdd = Equation_of_motion(params, th, phi, dth, dphi, u);

  % 3. Linearize the system at 0
  % You should end up with A (4x4), and b (4x1)
  x=[th phi dth dphi];
  x_dot=[dth dphi qdd(1) qdd(2)]';

  A=jacobian(x_dot,x);
  B=jacobian(x_dot,u);
  A=double(subs(A,[x u],zeros(1,5)));
  B=double(subs(B,[x u],zeros(1,5)));

  % 4. Check that (A,b) is  controllable
  % Number of uncontrollable states should return 0

  co=ctrb(A,B);
  unco=length(A)-rank(co);
 
% 5. Use LQR to get K as shown in the lecture
  K=lqr(A,B,eye(4),1,0);

  u=-K*X;
  
end

