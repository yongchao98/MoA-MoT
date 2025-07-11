#
# This script must be run with SageMath's Python interpreter, e.g., by running
# sage -python your_script_name.py
#
try:
    from sage.all import Knot, BraidGroup, var
except ImportError:
    print("This script requires the SageMath library.")
    print("Please install SageMath and run the script using 'sage -python <filename>'.")
    exit()

# === Step 1: Calculate the lower bound for Seifert circles of K1 ===
print("--- Analyzing K1 (the 10_74 knot) ---")

# Define the symbolic variables for the HOMFLY polynomial (a, z convention)
# This corresponds to the skein relation: a*P(L+) - a^{-1}*P(L-) = z*P(L0)
a, z = var('a, z')

# Get the knot K1 from the standard knot tables
K1 = Knot('10_74')

# Calculate the HOMFLY polynomial for K1
P1_homfly = K1.homfly_polynomial(normalization='az', variables=(a, z))

print(f"The HOMFLY polynomial for K1 (10_74) is: {P1_homfly}")

# To find the span in z, we find the min and max degrees of z in the polynomial
z_degrees = [m.degree(z) for m in P1_homfly.monomials()]
max_z_degree = max(z_degrees)
min_z_degree = min(z_degrees)

# The span is the difference between the max and min degrees
span_z = max_z_degree - min_z_degree

# The lower bound for the number of Seifert circles is span_z + 1
seifert_circle_lower_bound = span_z + 1

print(f"Maximum power of z in the polynomial: {max_z_degree}")
print(f"Minimum power of z in the polynomial: {min_z_degree}")
print(f"The z-span of the polynomial is {max_z_degree} - ({min_z_degree}) = {span_z}")
print(f"The lower bound of the minimum number of Seifert circles for K1 is {span_z} + 1 = {seifert_circle_lower_bound}.\n")


# === Step 2: Calculate the braid index of K2 ===
print("--- Analyzing K2 (closure of (sigma_1^-1)^3 * sigma_2^-1) ---")

# Define the 3-strand Braid Group
B3 = BraidGroup(3)
s1, s2 = B3.generators()

# Define the braid word for K2
braid_k2 = (s1**-1)**3 * s2**-1
print(f"The braid word for K2 is: {braid_k2}")

# Create the knot K2 by taking the closure of the braid
K2 = Knot(braid_k2)
# The knot K2 is actually the mirror of the 5_2 knot.
print(f"The knot K2 is identified as: {K2.name()}")

# Calculate the braid index of K2
braid_index_k2 = K2.braid_index()
print(f"The braid index of K2 is {braid_index_k2}.\n")

# === Step 3: Calculate the difference ===
print("--- Final Calculation ---")
difference = abs(braid_index_k2 - seifert_circle_lower_bound)

print(f"The difference between the braid index of K2 and the Seifert circle lower bound of K1 is:")
print(f"|{braid_index_k2} - {seifert_circle_lower_bound}| = {difference}")

<<<4>>>