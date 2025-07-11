from fractions import Fraction

# Define the rotation numbers for g and h based on the problem statement.
g_rot_num = 2
g_rot_den = 27
h_rot_num = 16
h_rot_den = 27

# Represent the rotation numbers as fractions for precision.
g_rot = Fraction(g_rot_num, g_rot_den)
h_rot = Fraction(h_rot_num, h_rot_den)

print("This program calculates the stable commutator length (scl) of g1*h2 in the free product G1*G2.")
print(f"The rotation number of g1 is rho(g1) = {g_rot.numerator}/{g_rot.denominator}.")
print(f"The rotation number of h2 is rho(h2) = {h_rot.numerator}/{h_rot.denominator}.")
print("-" * 30)

# The stable commutator length of an element f in G is taken to be (1/2) * |rho(f)|.
# Since the rotation numbers are positive, we omit the absolute value.
scl_g = Fraction(1, 2) * g_rot
scl_h = Fraction(1, 2) * h_rot

print(f"The scl of g1 in G1 is calculated as (1/2) * rho(g1) = {scl_g.numerator}/{scl_g.denominator}.")
print(f"The scl of h2 in G2 is calculated as (1/2) * rho(h2) = {scl_h.numerator}/{scl_h.denominator}.")
print("-" * 30)


# For an element w = g1*h2 in the free product G1*G2, the scl is given by
# scl(w) = max(scl(g1), scl(h2)).
result = max(scl_g, scl_h)

print("The final scl is the maximum of the individual scl values.")
# Output the final equation with all numbers, as requested.
print(f"scl(g1*h2) = max( (1/2) * ({g_rot.numerator}/{g_rot.denominator}), (1/2) * ({h_rot.numerator}/{h_rot.denominator}) )")
print(f"           = max( {scl_g.numerator}/{scl_g.denominator}, {scl_h.numerator}/{scl_h.denominator} )")
print(f"           = {result.numerator}/{result.denominator}")

# Final Answer
# print(f"The final computed stable commutator length is {result.numerator}/{result.denominator}.")
<<<8/27>>>