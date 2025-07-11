from fractions import Fraction

# Step 1: Define the rotation numbers for g and h.
rho_g = Fraction(2, 27)
rho_h = Fraction(16, 27)

# Step 2: Define the scl-norm of the rotation number quasimorphism rho.
# From the Bavard-Ghys theorem applied to G, scl(f) = (1/2)|rho(f)| for f in [G,G].
# ||rho||_scl = sup |rho(f)|/scl(f) = sup |rho(f)|/((1/2)|rho(f)|) = 2.
norm_rho_scl = 2

# Step 3: Compute the dual scl norms, scl_G^*(a) = |rho(a)| / ||rho||_scl.
scl_star_g = abs(rho_g) / norm_rho_scl
scl_star_h = abs(rho_h) / norm_rho_scl

# Step 4: Use the formula for scl in a free product: scl(g1*h2) = 1/2 * (scl_G^*(g1) + scl_G^*(h2)).
scl_final = Fraction(1, 2) * (scl_star_g + scl_star_h)

# Step 5: Print the calculation step-by-step.
# Note that scl_star_g is (|2/27|)/2 and scl_star_h is (|16/27|)/2.
# So we are computing 1/2 * ( (2/27)/2 + (16/27)/2 ).
term1_num = rho_g.numerator
term1_den = rho_g.denominator
term2_num = rho_h.numerator
term2_den = rho_h.denominator

print(f"scl(g1*h2) = 1/2 * (scl_G^*(g1) + scl_G^*(h2))")
print(f"           = 1/2 * ( (|rho(g)| / ||rho||) + (|rho(h)| / ||rho||) )")
print(f"           = 1/2 * ( |{term1_num}/{term1_den}|/{norm_rho_scl} + |{term2_num}/{term2_den}|/{norm_rho_scl} )")
print(f"           = 1/2 * ( {scl_star_g} + {scl_star_h} )")
print(f"           = 1/2 * ( {scl_star_g + scl_star_h} )")
print(f"           = {scl_final}")
