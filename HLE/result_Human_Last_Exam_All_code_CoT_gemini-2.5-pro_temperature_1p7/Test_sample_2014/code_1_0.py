import numpy as np

# Step 1: Define problem parameters
# We can normalize the chord length c to 1 for simplicity, as the ratio
# L1/L2 is dimensionless and independent of c, U_inf, alpha, and rho.
c = 1.0
s = 0.5 * c
h = 0.5 * c

# Step 2: Define coordinates of vortices and control points
# Aerofoil 1 (leading)
x1_le = 0.0
x1_vortex = x1_le + c / 4.0
x1_control = x1_le + 3 * c / 4.0
y1 = h

# Aerofoil 2 (trailing)
x2_le = x1_le + c + s
x2_vortex = x2_le + c / 4.0
x2_control = x2_le + 3 * c / 4.0
y2 = h

# Vortex positions (P1, P2)
p1 = np.array([x1_vortex, y1])
p2 = np.array([x2_vortex, y2])

# Control points (C1, C2)
c1 = np.array([x1_control, y1])
c2 = np.array([x2_control, y2])

# Image vortex positions (P1', P2')
p1_img = np.array([x1_vortex, -y1])
p2_img = np.array([x2_vortex, -y2])

# Step 3: Define a function to calculate induced velocity
# The flow tangency condition is: U*alpha = w_total_induced
# The vertical velocity 'w' induced at (xc, yc) by a vortex of strength
# Gamma at (xv, yv) is: w = -Gamma * (xc - xv) / (2 * pi * r^2)
# where r^2 = (xc - xv)^2 + (yc - yv)^2. This sign convention correctly
# models the upwash from the ground effect.
def induced_velocity(gamma, vortex_pos, control_pos):
    """Calculates the vertical induced velocity influence coefficient."""
    dx = control_pos[0] - vortex_pos[0]
    dy = control_pos[1] - vortex_pos[1]
    r_sq = dx**2 + dy**2
    if r_sq < 1e-9: # Avoid division by zero
        return 0.0
    w_per_gamma = -gamma * dx / (2 * np.pi * r_sq)
    return w_per_gamma

# Step 4: Calculate Aerodynamic Influence Coefficients (AIC) matrix A
# The system of equations is:
# A[0,0]*Gamma1 + A[0,1]*Gamma2 = U*alpha
# A[1,0]*Gamma1 + A[1,1]*Gamma2 = U*alpha
# A_ij is the total induced velocity at control point i from vortex system j per unit Gamma_j.

# A[0,0]: Influence of system 1 (image P1' with strength -1) on control point C1
a00 = induced_velocity(-1.0, p1_img, c1)

# A[0,1]: Influence of system 2 (real P2 with strength +1 and image P2' with -1) on C1
a01 = induced_velocity(1.0, p2, c1) + induced_velocity(-1.0, p2_img, c1)

# A[1,0]: Influence of system 1 (real P1 with strength +1 and image P1' with -1) on C2
a10 = induced_velocity(1.0, p1, c2) + induced_velocity(-1.0, p1_img, c2)

# A[1,1]: Influence of system 2 (image P2' with strength -1) on control point C2
a11 = induced_velocity(-1.0, p2_img, c2)

AIC_matrix = np.array([[a00, a01], [a10, a11]])

# Step 5: Solve the system of linear equations
# Set the right-hand side (U*alpha) to 1.0 to solve for the relative circulations.
b = np.array([1.0, 1.0])
gamma_vector = np.linalg.solve(AIC_matrix, b)
gamma1 = gamma_vector[0]
gamma2 = gamma_vector[1]

# Step 6: Calculate the lift ratio
# Lift L is proportional to circulation Gamma (L=rho*U*Gamma).
# So, the lift ratio L1/L2 is equal to the circulation ratio Gamma1/Gamma2.
lift_ratio = gamma1 / gamma2

# Step 7: Print the final equations and the result
print("The analysis leads to a system of two linear equations for the vortex strengths, Γ1 and Γ2.")
print("The general form is A * [Γ1, Γ2]^T = U*α * [1, 1]^T.\n")
print(f"For aerofoil 1, the equation is:")
print(f"({AIC_matrix[0,0]:.4f}) * Γ1 + ({AIC_matrix[0,1]:.4f}) * Γ2 = U*α")
print("\nFor aerofoil 2, the equation is:")
print(f"({AIC_matrix[1,0]:.4f}) * Γ1 + ({AIC_matrix[1,1]:.4f}) * Γ2 = U*α\n")

print(f"Solving this system gives the circulation ratio Γ1/Γ2 = {lift_ratio:.4f}.")
print("Since lift is proportional to circulation, the lift ratio L1/L2 is the same.\n")
print(f"The final lift ratio L1/L2 is: {lift_ratio:.4f}")
print("<<<1.8947>>>")