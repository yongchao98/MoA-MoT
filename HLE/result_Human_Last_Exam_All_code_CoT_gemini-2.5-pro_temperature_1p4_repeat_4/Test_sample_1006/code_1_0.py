# This script solves a topology problem by applying a uniqueness theorem.

# Step 1: Identify the space from the problem description.
# The problem describes a compact space X which is a compactification of the
# long ray R = [0, omega_1). The key property is that every bounded
# continuous function from R to the real numbers extends uniquely to X.

# Step 2: Recognize the defining property.
# The property of unique extension for all bounded continuous real-valued functions
# is the defining characteristic of the Stone-Čech compactification.
# The long ray R is a Tychonoff space, so its Stone-Čech compactification,
# denoted as βR, exists.

# Step 3: Apply the uniqueness theorem.
# The Stone-Čech compactification of a Tychonoff space is unique up to
# homeomorphism. This means any two spaces satisfying the given properties
# are homeomorphic to each other.

# Step 4: Determine the number of homeomorphism classes.
# Since all spaces satisfying the properties belong to the same homeomorphism class,
# there is only one such class.
number_of_classes = 1

# Final Equation: Number of classes = 1
print("The final equation is:")
print("Number of distinct homeomorphism classes = ", number_of_classes)
