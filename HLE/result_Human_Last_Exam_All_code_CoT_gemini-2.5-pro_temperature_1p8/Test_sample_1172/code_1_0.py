import sympy as sp

# Define symbolic variables
mu_0, I, h, R1, d = sp.symbols('mu_0 I h R_1 d', real=True, positive=True)

# Step 1: Define the source vector potential (dipole approximation)
# A_z_src = - (mu_0 * I * h * cos(theta)) / (2 * pi * r)
# The calculation relies on solving for the induced field.

# Step 2: Determine the induced potential coefficient K1
# The condition is that the total potential A_z(R1) is constant.
# A_z_tot = A_z_src + A_z_ind
# - (mu_0*I*h*cos(theta))/(2*pi*R1) + K1*R1*cos(theta) = const
# This implies the coefficient of cos(theta) is zero.
# K1*R1 - (mu_0*I*h)/(2*pi*R1) = 0
K1 = (mu_0 * I * h) / (2 * sp.pi * R1**2)

# Step 3: Find the induced magnetic field from the induced potential
# A_z_ind = K1 * r * cos(theta) = K1 * x
# B_ind is the curl of A_z_ind * z_hat.
# B_y_ind = -dA_z_ind/dx
B_y_ind = -K1

# Step 4: Calculate the flux of the induced field through the second circuit
# The second circuit has a width 'h'. The flux per unit length is B_y_ind * h.
Flux_ind_per_length = B_y_ind * h

# Step 5: Calculate the change in mutual inductance per unit length
# Delta_M = Flux_ind / I
Delta_M_per_length = Flux_ind_per_length / I

# Display the derivation and the final result
print("The change in mutual inductance per unit length, ΔM_l = M₂ - M₁, is due to the magnetic field induced by the concentrator shell.")
print("\n1. The ideal concentrator forces the total vector potential A_z to be constant on its inner surface (r=R₁).")
print(f"2. The vector potential from the source circuit (a dipole) is A_z_src ∝ (I*h*cos(θ))/r.")
print(f"3. The induced potential inside is A_z_ind = K₁*r*cos(θ).")
print(f"4. The boundary condition at r=R₁ allows solving for K₁:")
print(f"   K₁ = {sp.simplify(K1 / (I*h))} * I * h") # Show K1 in terms of variables
print(f"5. This induced potential corresponds to a uniform magnetic field in the y-direction:")
print(f"   B_y_ind = -K₁ = {sp.simplify(B_y_ind)}")
print(f"6. The flux from this induced field through the second circuit (width h) per unit length is:")
print(f"   Φ_ind/l = B_y_ind * h = {sp.simplify(Flux_ind_per_length)}")
print(f"\n7. Finally, the change in mutual inductance per unit length is this flux divided by the current I:")
print(f"   ΔM_l = (Φ_ind/l) / I = {sp.simplify(Delta_M_per_length)}")
print("\nThe final expression is:")
final_expression_latex = sp.latex(Delta_M_per_length)

# To fulfill the final output format instruction
# Print each term in the final equation.
# Expression is: - mu_0 * h^2 / (2 * pi * R1^2)
print("\nFinal Equation Breakdown:")
print(f"-1 * (mu_0 * h**2) / (2 * pi * R1**2)")

<<<print(f"- ( {mu_0} * {h}**2 ) / ( 2 * {sp.pi} * {R1}**2 )")>>>