import math

# Define the base function F and its derivatives.
# F(u, u_bar) = u * (1 - u)**2 * exp(-u_bar)
def F(u, u_bar):
    return u * (1 - u)**2 * math.exp(-u_bar)

# F_1 = dF/du
def F1(u, u_bar):
    return (1 - 4*u + 3*u**2) * math.exp(-u_bar)

# F_2 = dF/d(u_bar)
def F2(u, u_bar):
    return -u * (1 - u)**2 * math.exp(-u_bar)

# F_11 = d^2F/du^2
def F11(u, u_bar):
    return (-4 + 6*u) * math.exp(-u_bar)

# F_111 = d(F_11)/du
def F111(u, u_bar):
    return 6 * math.exp(-u_bar)

# d(F_11)/d(u_bar)
def F11_ubar(u, u_bar):
    return -(-4 + 6*u) * math.exp(-u_bar)

# We evaluate the expression at a point x, considering a shock at x+1.
# State at x: u0, u_bar0
# State at x+1: u1, u_bar1
# Configuration for maximum value:
u0 = 0.0
u_bar0 = 0.0
u1 = 1.0
u_bar1 = 1.0

# Evaluate derivatives at (u0, u_bar0)
F1_val = F1(u0, u_bar0)
F2_val = F2(u0, u_bar0)
F111_val = F111(u0, u_bar0)
F11_ubar_val = F11_ubar(u0, u_bar0)

# Evaluate F at (u0, u_bar0) and (u1, u_bar1)
F_u0_val = F(u0, u_bar0)
F_u1_val = F(u1, u_bar1)

# Calculate Lu and L_ubar
# Lu = -F2 * (u1 - u0)
Lu = -F2_val * (u1 - u0)

# L_ubar = F(u0, u_bar0) - F(u1, u_bar1) + F1(u0, u_bar0) * (u1 - u0)
L_ubar = F_u0_val - F_u1_val + F1_val * (u1 - u0)

# Calculate the final expression E = (d/dt + F1*d/dx)F11
# E = F111 * Lu + F11_ubar * L_ubar
E = F111_val * Lu + F11_ubar_val * L_ubar

# The final expression is an equation. We print all the components.
# E = F111(0,0) * [-F2(0,0)*(1-0)] + F11_ubar(0,0) * [F(0,0) - F(1,1) + F1(0,0)*(1-0)]
print(f"The maximum is obtained from the equation:")
print(f"E = {F111_val:.4f} * [-{F2_val:.4f} * ({u1} - {u0})] + {F11_ubar_val:.4f} * [{F_u0_val:.4f} - {F_u1_val:.4f} + {F1_val:.4f} * ({u1} - {u0})]")
print(f"E = {F111_val:.4f} * {Lu:.4f} + {F11_ubar_val:.4f} * {L_ubar:.4f}")
print(f"E = {F111_val * Lu:.4f} + {F11_ubar_val * L_ubar:.4f}")
print(f"Maximum value = {E}")
<<<4>>>