% Simulate the ideal case for the nondimensionalized ODE
function [t, v, vprime] = idealcase_sim(v0, vprime0, gamma, omega, tol)
  sys0 = [v0; vprime0]; % Vector of initial conditions
  while true
     [T, SYS] = ode45(@(t, sys) rhs(t, sys, gamma, omega), [0 (2*pi/omega)], sys0);
     dist = sqrt((SYS(end, 1) - sys0(1))^2 + (SYS(end, 2) - sys0(2))^2);
     if dist < tol
         break
     end
     sys0 = [SYS(end, 1); SYS(end, 2)];
  end
  t = T;
  v = SYS(:, 1);
  vprime = SYS(:, 2);
end
function sys_deriv = rhs(t, sys, gamma, omega) % sys(1) holds v and sys(2) holds v'
  sys_deriv = zeros(2, 1); % sys_deriv(1) holds v' and sys_deriv(2) holds v"
  sys_deriv(1) = sys(2);
  sys_deriv(2) = -2*gamma*sys(2) - sys(1) + sin(omega*t);
end
