% Simulate the nonideal case for the nondimensionalized ODE
function [t, v, u, R] = nonidealcase_sim(v0, u0, gamma, omega, epsilon, tol)
  sys0 = [v0; u0; 0];

  if isnan(tol)  % NaN signals: just integrate one period, don't preform a convergence check.
      [T, SYS] = ode45(@(t, sys) rhs(t, sys, gamma, omega, epsilon), [0 (2*pi/omega)], sys0, odeset('RelTol',1e-10,'AbsTol',1e-12));
  else
      max_iter = 2000;
      iter = 1;
      while iter < max_iter
          [T, SYS] = ode45(@(t, sys) rhs(t, sys, gamma, omega, epsilon), [0 (2*pi/omega)], sys0);
          dist = sqrt((SYS(end,1) - sys0(1))^2 + (SYS(end,2) - sys0(2))^2);
          if dist < tol, break; end
          sys0 = [SYS(end,1); SYS(end,2); 0];
          iter = iter + 1;
      end
  end

  t = T;
  v = SYS(:, 1);
  u = SYS(:, 2);
  A = SYS(end, 3);
  R = sqrt(abs(A)/(pi*omega));
end
function sys_deriv = rhs(t, sys, gamma, omega, epsilon) % sys(1) holds v and sys(2) holds u
  sys_deriv = zeros(3, 1); % sys_deriv(1) holds v' and sys_deriv(2) holds u'
  sys_deriv(1) = sys(2);
  sys_deriv(2) = -2*gamma*sys(2) - sys(1) - epsilon*(sys(1))^3 + sin(omega*t);
  sys_deriv(3) = 0.5*(sys(2)*sys_deriv(1) - sys(1)*sys_deriv(2));
end
