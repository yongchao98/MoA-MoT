# SageMath Code
# Define the braid group on 3 strands
B3 = BraidGroup(3)
s1, s2 = B3.gens()

# Define the braid for K2
beta = (s1**-1)**3 * s2**-1

# Create the knot K2 from the closure of the braid
K2 = Knot(beta)

# Get the braid index of K2
braid_index_K2 = K2.braid_index()

# Define the knot K1 = 10_74
K1 = Knot('10_74')

# Get the HOMFLY polynomial of K1 in v,z variables
P_K1 = K1.homfly_polynomial()

# Get the v-span of the polynomial
v_span = P_K1.degree(v) - P_K1.min_degree(v)

# Calculate the lower bound for the number of Seifert circles
lower_bound_s_K1 = v_span / 2 + 1

# Calculate the difference
difference = braid_index_K2 - lower_bound_s_K1

print(f"Braid for K2: {beta}")
print(f"Knot K2 is: {K2}")
print(f"Braid index of K2: {braid_index_K2}")
print(f"HOMFLY polynomial of K1 (10_74): {P_K1}")
print(f"v-span of P_K1: {v_span}")
print(f"Lower bound for Seifert circles of K1: {lower_bound_s_K1}")
print(f"Difference: {difference}")
