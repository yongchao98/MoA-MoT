try:
    import spherogram
except ImportError:
    print("Please install the 'spherogram' library by running: pip install spherogram")
    print("This library is required to compute knot polynomials.")
    # Exit gracefully if the library is not available
    exit()

# 1. Identify the knot.
# The knot 9_23 is the 23rd knot with 9 crossings in the standard knot tables.
knot_id = (9, 23)
knot = spherogram.Knot(*knot_id)

# 2. Compute the HOMFLY polynomial.
# Spherogram uses variables 'v' and 'z' for the HOMFLY polynomial P(v,z),
# which correspond to the standard 'a' and 'z' variables.
homfly_poly = knot.homfly_polynomial()

# The polynomial is stored in a dictionary-like structure where keys are tuples
# of exponents (exponent_v, exponent_z). We need the exponents for 'z'.
z_exponents = [exponent[1] for exponent in homfly_poly.poly.keys()]

# 3. Find the minimum and maximum powers of z.
j_min = min(z_exponents)
j_max = max(z_exponents)

# 4. Calculate the z-span of the polynomial.
span_z = j_max - j_min

# 5. Apply the formula to find the lower bound for the number of Seifert circles, s(K).
# The formula is s(K) >= (span_z / 2) + 1.
# Since span_z is always even for a knot, we can use integer division.
lower_bound = (span_z // 2) + 1

# 6. Print the results clearly, showing the full equation.
print(f"Finding the lower bound for the minimum number of Seifert circles, s(K), for the {knot} knot.")
print("-" * 70)
print(f"The HOMFLY polynomial P(v,z) for the {knot} knot is:")
print(f"P({knot}) = {homfly_poly}")
print()
print("From the polynomial, we extract the powers of the variable z.")
print(f"The minimum power of z is j_min = {j_min}")
print(f"The maximum power of z is j_max = {j_max}")
print()
print("The z-span of the polynomial is calculated as:")
print(f"span_z = j_max - j_min = {j_max} - {j_min} = {span_z}")
print()
print("A known theorem states that s(K) >= (span_z / 2) + 1.")
print("Substituting the values for the 9_23 knot, we get:")
print(f"s({knot}) >= ({span_z} / 2) + 1 = {span_z // 2} + 1 = {lower_bound}")
print("-" * 70)
print(f"The calculated lower bound for the minimum number of Seifert circles is {lower_bound}.")
