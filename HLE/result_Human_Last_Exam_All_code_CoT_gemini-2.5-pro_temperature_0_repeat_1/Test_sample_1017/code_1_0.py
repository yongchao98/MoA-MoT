from fractions import Fraction

# Define the rotation numbers for g and h
rot_g = Fraction(2, 27)
rot_h = Fraction(16, 27)

# The surface quasimorphism s(w) for w = g1*h2 is given by 1/2 * (rot(g1) - rot(h2))
# We take the absolute value as stable commutator length is non-negative.
s_val = Fraction(1, 2) * (rot_g - rot_h)
scl = abs(s_val)

# Print the calculation step-by-step
print("The stable commutator length is computed based on the surface quasimorphism associated with the rotation number.")
print("Let rot(g1) = 2/27 and rot(h2) = 16/27.")
print(f"scl(g1 h2) = |1/2 * (rot(g1) - rot(h2))|")
print(f"           = |1/2 * ({rot_g.numerator}/{rot_g.denominator} - {rot_h.numerator}/{rot_h.denominator})|")
print(f"           = |1/2 * ({(rot_g - rot_h).numerator}/{(rot_g - rot_h).denominator})|")
print(f"           = |{s_val.numerator}/{s_val.denominator}|")
print(f"           = {scl.numerator}/{scl.denominator}")
