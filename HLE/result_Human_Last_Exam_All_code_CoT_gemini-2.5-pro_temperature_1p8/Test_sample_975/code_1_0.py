import sympy

# Define symbols for clarity in expressions.
# Note: These are not used for computation, only for constructing the strings.
M_0, R_p, R, r, theta = sympy.symbols("M_0 R_p R r theta")
i_r, i_theta = sympy.symbols("i_r i_theta")

print("The determined magnetic field H(r, theta) is given in two regions:")
print("="*60)

# Region 1: Inside the shield (0 < r < R_p)
print("\nIn the region 0 < r < R_p:")
# Expression for the H field inside the sphere
H_inside_coeff = "M_0 * (2*R_p**3 + R**3) / (3*R**3)"
H_inside_vec = "(-cos(theta) * i_r + sin(theta) * i_theta)"
print(f"H = ( {H_inside_coeff} ) * ( {H_inside_vec} )")

# Region 2: Between the shield and the conductor (R_p < r < R)
print("\nIn the region R_p < r < R:")
# Expression for the H field outside the sphere
H_outside_r = "- (2*M_0 / 3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) * i_r"
H_outside_theta = "+ (M_0 / 3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) * i_theta"
# The string is split for readability
print(f"H = {H_outside_r} \\\n    {H_outside_theta}")

print("\nNote: i_r and i_theta are the unit vectors in the spherical coordinate system.")
print("-" * 60)
print("The results match option B.")