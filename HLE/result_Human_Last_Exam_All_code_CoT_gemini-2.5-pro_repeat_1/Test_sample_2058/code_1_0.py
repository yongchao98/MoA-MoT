import fractions

# The order of the group A_5
G_order = 60

# We calculate the total sum S = sum(2^(-c_rho(phi))) over all homomorphisms phi.
# This sum is broken down by the image of phi.

# Contribution from the trivial homomorphism (image {1})
S_trivial = 1
print(f"Contribution from trivial homomorphism: {S_trivial}")

# Contribution from homomorphisms with image C_2
# This is 15 * (3 * 2^(-4) + 4 * 2^(-6))
S_C2_val = 15 * (3 * (2**-4) + 4 * (2**-6))
S_C2_frac = fractions.Fraction(S_C2_val)
print(f"Contribution from C_2 images: 15 * (3/16 + 4/64) = {S_C2_frac.numerator}/{S_C2_frac.denominator}")

# Contribution from homomorphisms with image C_3
# This is 10 * (2^0 + 2^(-2))
S_C3_val = 10 * (1 + 1/4)
S_C3_frac = fractions.Fraction(S_C3_val)
print(f"Contribution from C_3 images: 10 * (1 + 1/4) = {S_C3_frac.numerator}/{S_C3_frac.denominator}")

# Contribution from homomorphisms with image V_4
# This is 30 * (3 * 2^(-6) + 6 * 2^(-8))
S_V4_val = 30 * (3 * (2**-6) + 6 * (2**-8))
S_V4_frac = fractions.Fraction(S_V4_val)
print(f"Contribution from V_4 images: 30 * (3/64 + 6/256) = {S_V4_frac.numerator}/{S_V4_frac.denominator}")

# Contribution from homomorphisms with image C_5
# This is 24 * 2^0
S_C5_val = 24 * 1
S_C5_frac = fractions.Fraction(S_C5_val)
print(f"Contribution from C_5 images: 24 * 1 = {S_C5_frac.numerator}")

# Total sum S
S_total_val = S_trivial + S_C2_val + S_C3_val + S_V4_val + S_C5_val
S_total_frac = fractions.Fraction(S_total_val)
print(f"Total sum S = {S_trivial} + {S_C2_frac} + {S_C3_frac} + {S_V4_frac} + {S_C5_frac} = {S_total_frac.numerator}/{S_total_frac.denominator}")

# Total mass M = S / |A_5|
M = S_total_frac / G_order
print(f"Total mass M = ({S_total_frac.numerator}/{S_total_frac.denominator}) / {G_order} = {M.numerator}/{M.denominator}")