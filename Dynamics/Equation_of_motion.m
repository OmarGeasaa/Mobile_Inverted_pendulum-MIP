function qdd = Equation_of_motion(params, th, phi, dth, dphi, u)
  
  % Provided params are
  % params.g: gravitational constant
  % params.mr: mass of the "rod"
  % params.ir: rotational inertia of the rod
  % params.d: distance of rod CoM from the wheel axis
  % params.r: wheel radius
  % Provided states are:
  % th: wheel angle (relative to body)
  % phi: body pitch
  % dth, dphi: time-derivatives of above
  % u: torque applied at the wheel
  
  q = [th;pi];
  x = params.r*(th+phi);
  
  % P is center of mass postion
  g = params.g;
  m = params.mr;
  i = params.ir;
  l = params.d;
  r = params.r;
  tau = u;
  
  A = zeros(2,2);
  A(1,1) = m*r^2;
  A(1,2) = m*r^2 + m*r*l*cos(phi);
  A(2,1) = m*r*l*cos(phi) + m*r^2;
  A(2,2) = 2*m*r*l*cos(phi) + i + m*l^2 + m*r^2;

  b = zeros(2,1);
  b(1,1) = m*r*l*dphi^2 * sin(phi) + tau;
  b(2,1) = m*l*sin(phi)*(g + r*dphi^2);
  
  qdd = A\b;
  % THE STUDENT WILL FILL THIS OUT
end