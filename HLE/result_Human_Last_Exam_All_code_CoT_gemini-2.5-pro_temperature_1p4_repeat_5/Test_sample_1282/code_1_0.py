import numpy as np

# The problem of finite-time blow-up can be analyzed by studying the linearized equation
# with a Gaussian initial condition, u_0(x) = exp(-a*x^2).
# The Fourier transform of the solution evolves according to:
# u_hat(k, t) = u_hat(k, 0) * exp(|k|^2 * (t + t^2/2))
# For u_0(x) = exp(-x^2), the Fourier transform is u_hat(k, 0) ~ exp(-k^2 / 4).
# So, u_hat(k, t) ~ exp(k^2 * (t + t^2/2 - 1/4)).
# Blow-up occurs when the coefficient of k^2 in the exponent becomes positive.
# The critical time T is found by setting the coefficient to zero.

# The equation to solve for the blow-up time T is:
# 0.5 * T^2 + T - 0.25 = 0

# We use the quadratic formula T = (-b +/- sqrt(b^2 - 4ac)) / (2a)
a = 0.5
b = 1.0
c = -0.25

# Calculate the discriminant
discriminant = b**2 - 4*a*c

print("We are solving the quadratic equation for the blow-up time T:")
print(f"{a}*T^2 + {b}*T + {c} = 0")
print("")

if discriminant < 0:
    print("No real solutions for the blow-up time.")
else:
    # Calculate the two possible solutions for T
    t1 = (-b + np.sqrt(discriminant)) / (2*a)
    t2 = (-b - np.sqrt(discriminant)) / (2*a)
    
    print(f"The mathematical solutions for T are {t1:.4f} and {t2:.4f}.")
    
    # The physically relevant solution must be positive
    if t1 > 0:
        blowup_time = t1
        print(f"The solution blows up at the finite positive time T = {blowup_time:.4f}")
    elif t2 > 0:
        blowup_time = t2
        print(f"The solution blows up at the finite positive time T = {blowup_time:.4f}")
    else:
        print("There is no positive blow-up time, so blow-up does not occur for this initial condition.")
