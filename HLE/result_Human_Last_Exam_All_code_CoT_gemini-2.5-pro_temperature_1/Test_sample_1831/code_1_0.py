# A program to calculate the number of equivalence classes of endomorphisms on a set of size 4.
# This is equivalent to counting the number of non-isomorphic functional graphs on 4 vertices.

# My method is to classify the endomorphisms by the number of their periodic points, k.
# n is the size of the set S.
n = 4

# We need the number of partitions of n, p(n), and the number of non-isomorphic rooted trees on n nodes, R(n).
# For the small values needed, we can use pre-calculated results.
# p(n): Number of integer partitions of n. p(2)=2, p(4)=5.
# R(n): Number of non-isomorphic rooted trees on n nodes. R(1)=1, R(2)=1, R(3)=2.

# Case 1: k = 4 periodic points.
# The function is a permutation. The number of classes is p(4).
num_classes_k4 = 5
print(f"Number of classes with {n} periodic points: {num_classes_k4}")

# Case 2: k = 3 periodic points.
# The permutation on the 3 periodic points can be a 3-cycle, a transposition, or 3 fixed points.
# The single transient point can be attached in different ways depending on the permutation's symmetry.
# Summing the possibilities gives 1 + 2 + 1 = 4 classes.
num_classes_k3 = 4
print(f"Number of classes with {n-1} periodic points: {num_classes_k3}")

# Case 3: k = 2 periodic points.
# The permutation on the 2 periodic points can be a 2-cycle or 2 fixed points.
# The 2 transient points form a forest. The number of ways to attach this forest gives the classes.
# For the 2-cycle, there are 3 classes. For the 2 fixed points, there are 3 classes. Total = 6.
num_classes_k2 = 6
print(f"Number of classes with {n-2} periodic points: {num_classes_k2}")

# Case 4: k = 1 periodic point.
# The function has one fixed point. The 3 transient points form a forest rooted at this point.
# The number of classes equals the number of non-isomorphic rooted forests on 3 nodes.
# This is R(3) + R(2)*R(1) + R(1)^3 = 2 + 1*1 + 1^3 = 4.
num_classes_k1 = 4
print(f"Number of classes with {n-3} periodic points: {num_classes_k1}")

# Total number of classes is the sum of the counts from each case.
total_classes = num_classes_k4 + num_classes_k3 + num_classes_k2 + num_classes_k1

print(f"\nTotal number of equivalence classes is the sum:")
print(f"{num_classes_k4} + {num_classes_k3} + {num_classes_k2} + {num_classes_k1} = {total_classes}")
