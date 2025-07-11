# A python program to calculate the number of equivalence classes.

# The total space X is a disjoint union of five component spaces.
# An equivalence relation x ~ y holds if there is a homeomorphism of X mapping x to y.

# As explained in the reasoning, any homeomorphism of X must map each of the five
# component spaces to itself, because no two are homeomorphic to each other.
# This means that equivalence classes are contained entirely within a single component space.
# We need to count the number of classes in each component and sum them up.

# For a space Y, the equivalence classes are the orbits of the points of Y
# under the action of the group of auto-homeomorphisms of Y.
# If this action is transitive, the space is called "homogeneous" and has 1 class.

# 1. The Torus
# The torus is homogeneous. For any two points x, y on the torus, there is a translation
# (which is a homeomorphism) that maps x to y.
num_classes_torus = 1

# 2. The Sphere
# The sphere is homogeneous. For any two points x, y on the sphere, there is a rotation
# (which is a homeomorphism) that maps x to y.
num_classes_sphere = 1

# 3. The Real Line
# The real line is homogeneous. The translation map f(z) = z + (y-x) is a homeomorphism
# mapping point x to point y.
num_classes_real_line = 1

# 4. A three-point discrete space
# In a discrete space, any bijection (permutation) of the points is a homeomorphism.
# The group of permutations on three points acts transitively.
num_classes_3_point_discrete = 1

# 5. A five-point discrete space
# Similarly, this space is homogeneous under the action of its permutation group.
num_classes_5_point_discrete = 1

# The total number of equivalence classes is the sum of the classes from each component.
total_classes = (num_classes_torus + num_classes_sphere + num_classes_real_line +
                 num_classes_3_point_discrete + num_classes_5_point_discrete)

# Print the final result, showing the numbers in the equation.
print(f"Number of equivalence classes on the torus: {num_classes_torus}")
print(f"Number of equivalence classes on the sphere: {num_classes_sphere}")
print(f"Number of equivalence classes on the real line: {num_classes_real_line}")
print(f"Number of equivalence classes on the three-point discrete space: {num_classes_3_point_discrete}")
print(f"Number of equivalence classes on the five-point discrete space: {num_classes_5_point_discrete}")
print("\nTotal number of equivalence classes is the sum:")
print(f"{num_classes_torus} + {num_classes_sphere} + {num_classes_real_line} + {num_classes_3_point_discrete} + {num_classes_5_point_discrete} = {total_classes}")

print("\n<<<5>>>")