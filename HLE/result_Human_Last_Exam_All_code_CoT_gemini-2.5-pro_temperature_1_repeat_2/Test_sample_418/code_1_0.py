# Each of the constituent spaces is homogeneous.
# A space is homogeneous if for any two points x and y, there exists
# an auto-homeomorphism of the space that sends x to y.
# This means that all points within a single homogeneous space belong to
# the same equivalence class.

# Number of equivalence classes in the torus
torus_classes = 1

# Number of equivalence classes in the sphere
sphere_classes = 1

# Number of equivalence classes in the real line
real_line_classes = 1

# Number of equivalence classes in a three-point discrete space
three_point_discrete_classes = 1

# Number of equivalence classes in a five-point discrete space
five_point_discrete_classes = 1

# The total space is a disjoint union of these five spaces.
# An auto-homeomorphism of the total space cannot map a point from one
# component space to another, as they are not homeomorphic.
# Therefore, the total number of equivalence classes is the sum of the
# number of classes in each component.
total_classes = (torus_classes +
                 sphere_classes +
                 real_line_classes +
                 three_point_discrete_classes +
                 five_point_discrete_classes)

# Print the final equation
print(f"{torus_classes} + {sphere_classes} + {real_line_classes} + {three_point_discrete_classes} + {five_point_discrete_classes} = {total_classes}")

# The final answer in the required format
print("<<<5>>>")