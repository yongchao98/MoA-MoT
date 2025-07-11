import sympy

# Define symbols for the parameters.
# g is the Yukawa coupling, epsilon is the dimensional regulator parameter.
# Since these will cancel out, we can represent the divergent term
# C = g^2 / (32 * pi^2 * epsilon) as a single symbolic constant.
C = sympy.Symbol('C')

print("Step 1: Determine counter-terms from the fermion self-energy Σ(p).")

# The divergent part of iΣ(p) in the MS-bar scheme is:
# iΣ_div = (g^2 / (16*pi^2*epsilon)) * ( (1/2)*p_slash + M_x )
# We can write this in terms of our constant C = g^2 / (32*pi^2*epsilon)
# iΣ_div = 2*C * ( (1/2)*p_slash + M_x ) = C*p_slash + 2*C*M_x
# Let C1 be the coefficient of p_slash and C2 be the constant term.
C1 = C
C2_over_Mx = 2*C
print(f"The divergent part of the fermion self-energy is iΣ_div = ({C1}) * p_slash + ({C2_over_Mx}) * M_x")

# The counter-term vertex is i*(δZ_x*p_slash - (δZ_x + δZ_m_x)*M_x)
# This must cancel -iΣ_div. So, the coefficients must match.
# δZ_x = -C1
# (δZ_x + δZ_m_x)*M_x = C2
delta_Zx = -C1
print(f"From the p_slash term, we get δZ_x = -C1 = {delta_Zx}")

# Now solve for delta_Z_mx
# δZ_m_x = C2/M_x - δZ_x
delta_Z_mx = C2_over_Mx - delta_Zx
print(f"From the mass term, (δZ_x + δZ_m_x) = C2/M_x = {C2_over_Mx}.")
print(f"Substituting δZ_x, we find δZ_m_x = {C2_over_Mx} - ({delta_Zx}) = {delta_Z_mx}")

print("\nStep 2: Determine δZ_g from the vertex correction Λ.")
# The divergent part of the vertex function Λ_div is -g^2/(32*pi^2*epsilon) = -C.
Lambda_div = -C
print(f"The divergent part of the vertex correction is Λ_div = {Lambda_div}")

# The relation between counter-terms is Λ_div = δZ_g + δZ_x (since δZ_φ = 0).
# δZ_g = Λ_div - δZ_x
delta_Zg = Lambda_div - delta_Zx
print(f"From the relation Λ_div = δZ_g + δZ_x, we find δZ_g = {Lambda_div} - ({delta_Zx}) = {delta_Zg}")

print("\nStep 3: Calculate the final ratio R.")
R = delta_Zx / (delta_Zg + delta_Z_mx)
print(f"R = δZ_x / (δZ_g + δZ_m_x)")
print(f"R = ({delta_Zx}) / (({delta_Zg}) + ({delta_Z_mx}))")
print(f"R = ({delta_Zx}) / ({delta_Zg + delta_Z_mx})")

# Final symbolic simplification
final_R = sympy.simplify(R)
print(f"\nThe final simplified value of the ratio is R = {final_R}.")

print(f"Let's print the numbers in the final equation: R = {delta_Zx.args[0] if delta_Zx.args else delta_Zx} / ({delta_Zg.args[0] if delta_Zg.args else delta_Zg} + {delta_Z_mx.args[0] if delta_Z_mx.args else delta_Z_mx}) = {final_R}")
