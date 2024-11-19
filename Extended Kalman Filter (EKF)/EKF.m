
function xhat = EKF(t, z)

  % vector of timestamps (of length T), and a 3xT matrix of observations, z.0
  n=length(t);
  xhat = zeros(2,length(t));
    % Student completes this
   
    Q=[5 0; 0 5];   % Q and R are tuning parameters
    R=[0.5 0 0; 0 0.5 0; 0 0 0.5];
    
    P=eye(2);    %prediction intial state
    for i=2:n
        dt=t(i)-t(i-1);
        A=[1 dt;0 1];
        xhat(:,i)=A*xhat(:,i-1);
        
        P=A*P*A.'+Q;
        phi=xhat(1,i);
        phi_dot=xhat(2,i);
        H = [cosd(phi)*(pi/180) -sind(phi)*(pi/180) 0; 0 0 1]'; 
        h = [sind(phi);cosd(phi);phi_dot];
        
        K = P*H.'/(H*P*H.'+R);
        xhat(:,i) =  xhat(:,i)+ K * (z(:,i) - h);
        P = (eye(2) - K*H) * P; 
        
    end

end
