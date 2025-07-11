# The problem asks for the smallest number of topologically distinct compactifications
# of the ray with a remainder X, where X is a nondegenerate, locally-connected,
# compact metric space.

# According to a theorem in topology, the number of such compactifications is
# equal to the number of non-empty, closed, connected subsets of X.

# Our goal is to find the minimum possible value for this number.

# Step 1: Establish a lower bound.
# The space X is 'nondegenerate', which means it must have at least two points.
# Let's call two distinct points in X, p1 and p2.
# In any topological space, a singleton set {p} is always closed and connected.
# Therefore, {p1} and {p2} are two distinct non-empty, closed, connected subsets of X.
# This means that any valid space X must have at least 2 such subsets.
# So, the minimum number is at least 2.

# Step 2: Find a space X that achieves this minimum.
# Consider the discrete space X with two points, e.g., X = {0, 1}.
# Let's check if this space satisfies the given conditions:
# - It is nondegenerate (has 2 points).
# - It is a compact metric space (being finite).
# - It is locally connected (the singleton neighborhoods are connected).

# Now, let's count the non-empty, closed, connected subsets of X = {0, 1}.
# - The non-empty subsets are {0}, {1}, and {0, 1}.
# - In a discrete space, all subsets are closed.
# - The connected subsets are only the singletons: {0} and {1}. The set {0, 1} is not connected.
# So, X = {0, 1} has exactly two non-empty, closed, connected subsets.

# Step 3: Conclusion.
# We have shown that the minimum number must be at least 2, and we have found a
# space for which the number is exactly 2.
# Therefore, the smallest possible number is 2.

smallest_number_of_compactifications = 2
print(smallest_number_of_compactifications)