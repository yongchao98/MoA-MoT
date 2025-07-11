import spherogram
import sympy

# Step 1: Define the knot and compute its HOMFLY polynomial.
# The knot is specified by its Alexander-Briggs notation '9_23'.
try:
    knot = spherogram.Knot('9_23')
except ImportError as e:
    print(f"Error: {e}")
    print("Please ensure you have the necessary data files for spherogram. They are usually included with the installation.")
    exit()

homfly_poly = knot.homfly_polynomial()

print(f"The HOMFLY polynomial P(a, z) for the 9_23 knot is:")
print(f"P(a, z) = {homfly_poly}\n")

# Step 2: Find the maximum and minimum degrees of the variable z.
# We treat the expression as a polynomial in z and find its term's powers.
poly_in_z = sympy.poly(homfly_poly, sympy.Symbol('z'))
z_powers = [term[0][0] for term in poly_in_z.terms()]

max_deg_z = max(z_powers)
min_deg_z = min(z_powers)

print(f"The maximum degree of z in the polynomial is: {max_deg_z}")
print(f"The minimum degree of z in the polynomial is: {min_deg_z}\n")

# Step 3: Calculate the span of z.
span_z = max_deg_z - min_deg_z
print(f"The span of z is calculated as: span_z = max_deg_z - min_deg_z")
print(f"span_z = {max_deg_z} - ({min_deg_z}) = {span_z}\n")

# Step 4: Use the Morton-Franks-Williams inequality to find the lower bound for s(K).
# The inequality is s(K) >= span_z / 2 + 1.
lower_bound = (span_z / 2) + 1

print("From the Morton-Franks-Williams inequality, we find the lower bound for the minimum number of Seifert circles, s(K):")
print("s(K) >= span_z / 2 + 1")
print(f"s(K) >= {span_z} / 2 + 1")
print(f"s(K) >= {span_z / 2} + 1")
print(f"s(K) >= {int(lower_bound)}")
print(f"\nTherefore, a lower bound for the minimum number of Seifert circles is {int(lower_bound)}.")
