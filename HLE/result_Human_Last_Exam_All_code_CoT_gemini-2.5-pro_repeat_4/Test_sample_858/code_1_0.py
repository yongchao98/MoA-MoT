# Step 1: State the problem and definitions.
# The problem asks for the smallest possible cardinality of the set of non-block points
# in an aposyndetic continuum X.

# Definitions:
# - Continuum (X): A compact, connected, Hausdorff space. By convention, continua are non-empty.
# - Aposyndetic (X): For every two distinct points x, y in X, there exists a subcontinuum K
#   such that x is in the interior of K (Int(K)) and K is a subset of X \ {y}.
# - Non-block point (p): A point p in X such that X \ {p} contains a dense continuum-connected subset.

# Step 2: Use a key theorem to simplify the problem.
# A fundamental theorem in continuum theory by F. B. Jones (1951) states:
# If a continuum X is aposyndetic, then every point of X is a non-block point.
#
# This is because for an aposyndetic continuum X, the subspace X \ {p} is continuum-connected for any p in X.
# A continuum-connected space is, by definition, its own dense continuum-connected subset.
# Therefore, the set of non-block points in an aposyndetic continuum is the entire space X.

# Step 3: Rephrase the question.
# The problem is now to find the smallest possible cardinality of an aposyndetic continuum X.
# Let N be the set of non-block points. The theorem implies N = X.
# We need to find the minimum possible value for |N| = |X|.

# Step 4: Find the minimum cardinality for an aposyndetic continuum X.
# A continuum must be a non-empty space, so its cardinality must be at least 1.
# We need to check if a space with cardinality 1 can be an aposyndetic continuum.
#
# Let's consider the space X = {p}, which contains a single point.
# - Is X a continuum? Yes. A single-point space is non-empty, compact, connected, and Hausdorff.
# - Is X aposyndetic? The definition of aposyndesis is a condition "for every two distinct points x, y in X...".
#   In the space X = {p}, there are no pairs of distinct points.
#   Therefore, the condition is vacuously true.
#   So, X = {p} is an aposyndetic continuum.

# Step 5: Conclude the result.
# We have found an aposyndetic continuum, X = {p}, with cardinality 1.
# Since the cardinality of any continuum must be at least 1, the minimum possible cardinality is 1.

smallest_cardinality = 1

# The problem asks to output the number in the final equation.
# Here, the final answer is simply the number we found.
print("The smallest possible cardinality of the set of non-block points is:")
print(smallest_cardinality)