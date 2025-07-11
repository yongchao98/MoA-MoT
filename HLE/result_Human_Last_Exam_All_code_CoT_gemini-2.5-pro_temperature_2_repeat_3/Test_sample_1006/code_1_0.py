# The user wants to find the number of distinct homeomorphism classes for a
# compact topological space X with specific properties related to the long ray R.

# Step 1: Analyze the properties of the space X.
# The problem states that X is a compact topological space such that:
# (1) X contains a dense copy of the long ray R = [0, omega_1). This means X is a compactification of R.
# (2) Every bounded continuous function f from R to the real numbers R extends to a unique continuous function from X to R.

# Step 2: Identify the space X using its properties.
# The two properties given are the defining characteristics of the Stone-Čech compactification of R, denoted as beta R.
# The long ray R is a Tychonoff space, so its Stone-Čech compactification exists.
# The universal property of the Stone-Čech compactification beta R is that for any continuous map g from R to a compact
# Hausdorff space K, there exists a unique continuous extension g_bar from beta R to K.
# Property (2) is a special case of this universal property where K is a closed, bounded interval in R (which is compact and Hausdorff).
# This special case is sufficient to uniquely characterize the Stone-Čech compactification.
# Therefore, any space X that satisfies the given conditions must be homeomorphic to beta R.

# Step 3: Determine the number of homeomorphism classes.
# The Stone-Čech compactification is unique up to homeomorphism. This means if we have two spaces, X1 and X2, that
# both satisfy the given properties, then both X1 and X2 must be homeomorphic to beta R.
# If X1 is homeomorphic to beta R and X2 is homeomorphic to beta R, then X1 must be homeomorphic to X2.
# This implies that all such spaces X belong to a single, unique homeomorphism class.

# Step 4: State the final answer.
# Since all spaces satisfying the conditions belong to the same homeomorphism class, there is only one such class.
# The number is 1.

final_answer = 1
print(f"The number of distinct homeomorphism classes is: {final_answer}")