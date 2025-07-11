# The reasoning for the number of homeomorphism classes has been explained above.
# This code block simply prints the final calculation as requested.

# The first class is for Peano continua, which must be homeomorphic to the arc [0, 1].
num_classes_peano = 1

# The second class is for non-Peano continua, exemplified by the topologist's sine curve.
# Advanced results show this is the only other possibility.
num_classes_non_peano = 1

# The total number of distinct homeomorphism classes is the sum of these cases.
total_classes = num_classes_peano + num_classes_non_peano

# Print the final equation, showing each number.
print(f"The number of classes for Peano continua is {num_classes_peano}.")
print(f"The number of classes for non-Peano continua is {num_classes_non_peano}.")
print(f"The total number of distinct homeomorphism classes is:")
print(f"{num_classes_peano} + {num_classes_non_peano} = {total_classes}")
