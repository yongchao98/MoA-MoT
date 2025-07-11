import math

# Step 1: Analyze the given properties.
# Let X be the continuum.
# Let E be the set of endpoints and I be the set of non-endpoints (the "interior").
# Let n be the number of endpoints.

# Property (1): 1 < n < infinity
# Property (2): X has exactly 2 orbits under its auto-homeomorphism group.

# Step 2: Deduce the structure of X.
# From Property (2), the two orbits must be E and I.
# This implies that the interior I is a homogeneous space (all its points are topologically equivalent).
# For a graph-like continuum, the homogeneity of I means there can be no branch points (vertices of degree > 2).
# This restricts X to be either a simple path (arc) or a simple cycle (circle).

# Step 3: Evaluate the candidates.
# A circle has n = 0 endpoints, which violates Property (1).
# An arc has n = 2 endpoints. This satisfies Property (1) because 2 > 1 and 2 is finite.

# Step 4: Verify the arc against Property (2).
# An arc (e.g., the interval [0,1]) has two orbits:
# Orbit 1: The set of endpoints {0, 1}.
# Orbit 2: The set of interior points (0, 1).
# The number of orbits is 2, which satisfies Property (2).

# Step 5: Conclude the number of distinct topological types.
# Any continuum with the given properties must be topologically equivalent to a simple arc.
# All simple arcs are homeomorphic to each other.
# Therefore, there is only one such topological type.

number_of_distinct_continua = 1

print("Deduction process:")
print("1. The two properties imply the continuum is partitioned into two orbits: the set of endpoints (E) and the set of non-endpoints (I).")
print("2. The set of non-endpoints (I) must be a homogeneous space.")
print("3. This structural requirement forces the continuum to have no 'branch points', meaning it must be topologically a simple path (an arc) or a simple cycle (a circle).")
print("4. A circle has 0 endpoints, which fails Property (1) (number of endpoints > 1).")
print("5. A simple arc has 2 endpoints, which satisfies Property (1).")
print("6. A simple arc has 2 orbits (the two endpoints, and the interior), which satisfies Property (2).")
print("7. Conclusion: The continuum must be a simple arc.")
print("\nSince all simple arcs are topologically equivalent, there is only one such topologically distinct continuum.")

final_equation = f"Number of topologically distinct continua = {number_of_distinct_continua}"
print(final_equation)

<<<1>>>