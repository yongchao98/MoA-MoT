# The problem of classifying homogeneous planar continua is solved by
# dividing them into two types based on whether they separate the plane.

# 1. According to a theorem by R. H. Bing, the number of homogeneous
#    planar continua that DO NOT separate the plane is 2.
#    (These are the circle and the pseudo-arc).
num_non_separating_classes = 2

# 2. According to a theorem by J. Prajs and J. T. Rogers, Jr., the number of
#    homogeneous planar continua that DO separate the plane is 3.
num_separating_classes = 3

# The total number of classes is the sum of these two counts.
total_classes = num_non_separating_classes + num_separating_classes

# Print the final equation and the answer.
print(f"The total number of homeomorphism classes of homogeneous planar continua is the sum of the non-separating and separating classes.")
print(f"{num_non_separating_classes} + {num_separating_classes} = {total_classes}")