import itertools

# The first 12 vectors (L1) are of the form e_i +/- e_j
# Let's represent e_i as a tuple with 1 at index i-1 and 0 otherwise.
e = [
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1)
]

# Generate the 12 vectors of L1
L1 = []
for i in range(4):
    for j in range(i + 1, 4):
        # e_i + e_j
        v_plus = tuple(e[i][k] + e[j][k] for k in range(4))
        # e_i - e_j
        v_minus = tuple(e[i][k] - e[j][k] for k in range(4))
        L1.append(v_plus)
        L1.append(v_minus)

# Print L1 to verify (optional)
# print("Vectors in L1:")
# for v in L1:
#    print(v)

# Count the number of linearly dependent triples.
# For vectors with integer components, linear dependence over H can be tested over Q.
# A set {v1, v2, v3} is dependent if the matrix with these rows has rank < 3.
# Or if one is a linear combination of the others: v3 = a*v1 + b*v2.
# A simpler way for these specific vectors: check for dependencies like v1 +/- v2 +/- v3 = 0
dependent_triples = set()
for combo in itertools.combinations(L1, 3):
    v1, v2, v3 = combo
    # Check for all 4 sign combinations for c1*v1 + c2*v2 + c3*v3 = 0
    # where ci are +/- 1. We only need to check one set of signs like (1,1,1)
    # as other relations are equivalent, e.g., v1+v2-v3=0 is same as v3=v1+v2.
    for signs in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]:
        is_dependent = True
        for i in range(4):
            if signs[0]*v1[i] + signs[1]*v2[i] + signs[2]*v3[i] != 0:
                is_dependent = False
                break
        if is_dependent:
            # Use frozenset to store the unique combination of vectors
            dependent_triples.add(frozenset(combo))
            break # Found a dependency, move to next combination

num_dependent_triples = len(dependent_triples)

print("The final equation is composed of the following numbers.")
print(f"Number of distinct indices to choose for a triangle: {4}")
print(f"Number of ways to choose signs for dependence relation: {4}")
print(f"Total number of dependent triples: {4} * {4} = {num_dependent_triples}")
print(f"Dimension of H^9(M,Q) = {num_dependent_triples}")