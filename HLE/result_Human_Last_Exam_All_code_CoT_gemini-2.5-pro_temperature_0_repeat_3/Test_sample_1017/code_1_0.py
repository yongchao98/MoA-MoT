from fractions import Fraction

# Define the translation values for g and h
g_val = Fraction(2, 27)
h_val = Fraction(16, 27)

# The stable commutator length (scl) is given by the formula |rot(g) - rot(h)| / 4
# where rot(g) and rot(h) are the translation values.

# Step 1: Calculate the difference in rotation numbers
diff = g_val - h_val
abs_diff = abs(diff)

# Step 2: Calculate the scl
scl = abs_diff / 4

# Print the steps of the calculation
print(f"The rotation number of g is rot(g) = {g_val.numerator}/{g_val.denominator}")
print(f"The rotation number of h is rot(h) = {h_val.numerator}/{h_val.denominator}")
print(f"The formula for the stable commutator length is scl = |rot(g) - rot(h)| / 4")
print(f"scl = |{g_val.numerator}/{g_val.denominator} - {h_val.numerator}/{h_val.denominator}| / 4")
print(f"scl = |{diff.numerator}/{diff.denominator}| / 4")
print(f"scl = ({abs_diff.numerator}/{abs_diff.denominator}) / 4")
print(f"scl = {scl.numerator}/{scl.denominator}")
