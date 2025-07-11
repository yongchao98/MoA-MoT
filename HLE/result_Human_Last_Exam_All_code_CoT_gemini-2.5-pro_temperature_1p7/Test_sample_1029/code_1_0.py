import sympy

# Set up sympy for pretty printing
sympy.init_printing(use_unicode=True)

# Define symbolic variables used in the physics problem
# rho: volume charge density
# E: magnitude of the external electric field
# v: speed of the rod
# R: radius of the rod
# r: radial distance from the axis of the rod
# mu_0: permeability of free space
# epsilon_0: permittivity of free space
rho, E, v, R, r, mu_0, epsilon_0 = sympy.symbols('rho E v R r mu_0 epsilon_0', real=True, positive=True)

# Define string representations for the basis vectors in a cylindrical coordinate system
r_hat_str = "r_hat"  # Radial unit vector
k_hat_str = "z_hat"  # Axial unit vector (along the rod's axis)

# --- Calculation for the region INSIDE the rod (radial distance r < R) ---

# Step 1: Define the electric field E_total = E_rod + E_ext
# E_rod (from Gauss's Law) has a radial component. E_ext is axial.
E_r_inside = (rho * r) / (2 * epsilon_0)
E_z_inside = E

# Step 2: Define the magnetic field B (from Ampere's Law)
# The moving charge creates a current, resulting in an azimuthal magnetic field.
B_phi_inside = (mu_0 * rho * v * r) / 2

# Step 3: Compute the Poynting vector S = (1/mu_0) * (E x B)
# The cross product in cylindrical coordinates yields:
# E x B = (E_r*B_phi) * (r_hat x phi_hat) + (E_z*B_phi) * (z_hat x phi_hat)
#       = (E_r*B_phi) * z_hat - (E_z*B_phi) * r_hat
# We compute the coefficients for the radial (r_hat) and axial (z_hat) components.
S_r_inside_comp = - (1 / mu_0) * (E_z_inside * B_phi_inside)
S_z_inside_comp = (1 / mu_0) * (E_r_inside * B_phi_inside)


# --- Calculation for the region OUTSIDE the rod (radial distance r > R) ---

# Step 1: Define the electric field E_total
E_r_outside = (rho * R**2) / (2 * epsilon_0 * r)
E_z_outside = E

# Step 2: Define the magnetic field B
B_phi_outside = (mu_0 * rho * v * R**2) / (2 * r)

# Step 3: Compute the Poynting vector components using the same logic
S_r_outside_comp = - (1 / mu_0) * (E_z_outside * B_phi_outside)
S_z_outside_comp = (1 / mu_0) * (E_r_outside * B_phi_outside)


# --- Format and Print the Final Results ---

# Use sympy.simplify() to get the simplest form of the expressions
s_inside_radial_final = sympy.simplify(S_r_inside_comp)
s_inside_axial_final = sympy.simplify(S_z_inside_comp)
s_outside_radial_final = sympy.simplify(S_r_outside_comp)
s_outside_axial_final = sympy.simplify(S_z_outside_comp)

# Output the final equations and their components clearly
print("The Poynting vector S is computed for two regions:")
print("=" * 60)

print("\n1. Inside the rod (for r < R):\n")
print(f"S_inside = ({s_inside_radial_final}) * {r_hat_str}   +   ({s_inside_axial_final}) * {k_hat_str}")
print("\n   where each component's coefficient is:")
print(f"   Radial ({r_hat_str}) coefficient: {s_inside_radial_final}")
print(f"   Axial ({k_hat_str}) coefficient:   {s_inside_axial_final}")

print("\n" + "=" * 60)

print("\n2. Outside the rod (for r > R):\n")
print(f"S_outside = ({s_outside_radial_final}) * {r_hat_str}   +   ({s_outside_axial_final}) * {k_hat_str}")
print("\n   where each component's coefficient is:")
print(f"   Radial ({r_hat_str}) coefficient: {s_outside_radial_final}")
print(f"   Axial ({k_hat_str}) coefficient:   {s_outside_axial_final}")
