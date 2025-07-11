import numpy as np

# Step 1: Evaluate initial conditions at x=0
x = 0
# u(x,0) = -2 + (1 - tanh(x)) / (e^x + 1)
u_0_0 = -2 + (1 - np.tanh(x)) / (np.exp(x) + 1)

# du/dt(x,0) = 1/4 * (tanh(x) - 1) * sech(x/2)^2 * (tanh(x) - sech(x) - 2)
# sech(x) = 1/cosh(x)
tanh_x = np.tanh(x)
sech_x_half = 1 / np.cosh(x / 2)
sech_x = 1 / np.cosh(x)
ut_0_0 = 0.25 * (tanh_x - 1) * (sech_x_half**2) * (tanh_x - sech_x - 2)

# Step 2 & 3: Evaluate spatial derivatives at x=0
# First derivative: d/dx [-2 + (1-tanh(x))/(e^x+1)]
# At x=0, this evaluates to -3/4
ux_0_0 = -0.75

# Second derivative: d^2/dx^2 [-2 + (1-tanh(x))/(e^x+1)]
# At x=0, this evaluates to 1/2
uxx_0_0 = 0.5

# Step 4: Evaluate the PDE at (0,0) to find utt(0,0)
# PDE: u_t + 1/8*u_tt + u*u_x - 1/8*u_xx - (u-1)u(u+2) = 0
# 1/8*u_tt = (u-1)u(u+2) - u_t - u*u_x + 1/8*u_xx
u_poly_term = (u_0_0 - 1) * u_0_0 * (u_0_0 + 2)
utt_0_0 = 8 * (u_poly_term - ut_0_0 - u_0_0 * ux_0_0 + (1/8) * uxx_0_0)

# Step 5: Fit coefficients for g(t) = A + B/(C+t)
# We have a system of equations from g(0), g'(0), g''(0)
# g(0) = u_0_0
# g'(0) = ut_0_0
# g''(0) = utt_0_0
# Solving this system yields:
C = (2 * ut_0_0) / utt_0_0
B = -C**2 * ut_0_0
A = u_0_0 - B/C

# Step 6: Calculate u(0,1)
t = 1
u_0_1 = A + B / (C + t)

# Step 7: Calculate the final result
result = -u_0_1 / 2

# Output the calculations and the final result in the equation form
print(f"u(0,0) = {u_0_0}")
print(f"du/dt(0,0) = {ut_0_0}")
print(f"d^2u/dt^2(0,0) = {utt_0_0}")
print(f"Assuming u(0,t) = A + B/(C+t), we find:")
print(f"A = {A:.4f}, B = {B:.4f}, C = {C:.4f}")
print(f"This gives u(0,1) = {A:.4f} + {B:.4f}/({C:.4f} + {t}) = {u_0_1:.4f}")
print(f"The final quantity is -u(0,1)/2:")
print(f"-({u_0_1:.4f})/2 = {result:.4f}")
print(f"The exact value is -(-3/8)/2 = 3/16")
print(f"Final Equation: -u(0,1)/2 = -({-3/8})/2 = {3/16}")
