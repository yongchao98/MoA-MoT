# The number of isomorphism classes of automorphism groups for a compact
# Riemann surface of a given genus is a result from the classification
# of finite group actions on surfaces. These values are well-established in
# mathematical literature.

# For a Riemann surface of genus g=2, the number of isomorphism classes
# of possible automorphism groups is 12.
num_groups_g2 = 12

# For a Riemann surface of genus g=3, the number of isomorphism classes
# of possible automorphism groups is 36.
num_groups_g3 = 36

# For a Riemann surface of genus g=4, the number of isomorphism classes
# of possible automorphism groups is 23.
num_groups_g4 = 23

# The problem requires the answer in the format [g2, g3, g4].
result = [num_groups_g2, num_groups_g3, num_groups_g4]

# We print the final list, which contains each of the required numbers.
print(result)