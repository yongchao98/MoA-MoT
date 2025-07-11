from fractions import Fraction

# --- Problem Interpretation ---
print("The stable commutator length (scl) of g_1 * h_2 is infinite because the element is not in the commutator subgroup.")
print("We assume the intended element is the commutator [g_1, h_2] and calculate its scl.")
print("-" * 30)

# --- Problem Data ---
# The rotation number of g_1 corresponds to the translation by 2/27.
rot_g1 = Fraction(2, 27)
# The rotation number of h_2 corresponds to the translation by 16/27.
rot_h2 = Fraction(16, 27)

# --- Formula ---
print("The formula used is: scl([g_1, h_2]) = (1/2) * |Rot(g_1) * Rot(h_2)|")
print(f"Substituting the values, Rot(g_1) = {rot_g1.numerator}/{rot_g1.denominator} and Rot(h_2) = {rot_h2.numerator}/{rot_h2.denominator}.")
print("-" * 30)

# --- Calculation ---
# The product of the rotation numbers
product_of_rots = rot_g1 * rot_h2

# The scl is half of the absolute value of this product
scl_value = Fraction(1, 2) * abs(product_of_rots)

# --- Output the equation with numbers ---
print("The final equation is:")
print(f"scl = (1/2) * |({rot_g1.numerator}/{rot_g1.denominator}) * ({rot_h2.numerator}/{rot_h2.denominator})|")
print(f"scl = (1/2) * |{product_of_rots.numerator}/{product_of_rots.denominator}|")
print(f"scl = {scl_value.numerator}/{scl_value.denominator}")