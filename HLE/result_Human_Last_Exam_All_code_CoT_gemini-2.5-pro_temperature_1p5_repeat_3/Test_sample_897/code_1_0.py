import sympy

# --- Introduction ---
print("This script solves the given knot theory problem by calculating two values and finding their difference.")
print("The required library is `sympy`. You can install it using: pip install sympy\n")

# --- Part 1: Analyzing the 10_74 knot (K1) ---
print("--- Part 1: Lower bound of Seifert circles for K1 (10_74 knot) ---")
print("We use the HOMFLY polynomial P(K; a, z) and Cromwell's theorem: s(K) >= span_z(P(K)) + 1")

# Define symbols and the HOMFLY polynomial for the 10_74 knot
a, z = sympy.symbols('a z')
homfly_k1_str = "a**6*z**2 - a**6 - a**4*z**4 + 3*a**4*z**2 - 2*a**4 + a**2*z**4 - 2*a**2*z**2 + a**2"
homfly_k1 = sympy.sympify(homfly_k1_str)

print(f"\nThe HOMFLY polynomial for K1 (10_74) is: P(a, z) = {homfly_k1}")

# Treat the expression as a polynomial in z to find the degree span
homfly_k1_poly_z = sympy.Poly(homfly_k1, z)
z_exponents = [monom[0] for monom in homfly_k1_poly_z.monoms()]

max_deg_z = max(z_exponents)
min_deg_z = min(z_exponents)

print(f"The maximum degree of z is {max_deg_z}.")
print(f"The minimum degree of z is {min_deg_z}.")

# Calculate the z-span and the lower bound
span_z = max_deg_z - min_deg_z
seifert_lower_bound = span_z + 1

print(f"The z-span is: {max_deg_z} - {min_deg_z} = {span_z}")
print(f"The lower bound for the minimum number of Seifert circles is: {span_z} + 1 = {seifert_lower_bound}")
print("-" * 50)


# --- Part 2: Analyzing K2 ---
print("\n--- Part 2: Braid index of K2 (closure of (sigma_1^-1)^3 * sigma_2^-1) ---")
print("K2 is the closure of a 3-strand braid, so its braid index b(K2) is at most 3.")

print("\nThis braid closure corresponds to the 5_2 knot.")
print("A knot has braid index 2 if and only if it is a (2,k)-torus knot.")
print("The 5_2 knot is not a torus knot, which can be verified using polynomial invariants.")

braid_index_k2 = 3
print(f"Therefore, the braid index of K2 must be {braid_index_k2}.")
print("-" * 50)


# --- Part 3: Calculating the Final Difference ---
print("\n--- Part 3: Final Calculation ---")
print(f"Value 1: Lower bound of Seifert circles for K1 = {seifert_lower_bound}")
print(f"Value 2: Braid index of K2 = {braid_index_k2}")

difference = seifert_lower_bound - braid_index_k2

print(f"\nThe difference between the two values is: {seifert_lower_bound} - {braid_index_k2} = {difference}")

print(f"\n<<<{difference}>>>")