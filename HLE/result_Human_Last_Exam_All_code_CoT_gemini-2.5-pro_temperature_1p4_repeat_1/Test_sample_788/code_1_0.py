# The number of equivalence classes is determined by the number of independent invariants.
# As explained in the reasoning, there are 3 independent binary invariants.
# The number of classes is 2 to the power of the number of independent invariants.

num_independent_invariants = 3
num_classes = 2**num_independent_invariants

print("The number of equivalence classes is 2 raised to the power of the number of independent binary invariants.")
print("Number of independent invariants = " + str(num_independent_invariants))
print(f"2^{num_independent_invariants} = {num_classes}")
