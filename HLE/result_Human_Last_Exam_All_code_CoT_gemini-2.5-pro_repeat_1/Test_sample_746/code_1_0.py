# The user wants to find the maximum number of dispersion points in a
# compact connected metric space. This is a problem from point-set topology.
# The solution involves a proof by contradiction.

# Step 1: Define the key concepts.
# A 'continuum' is a compact connected metric space.
# A point 'x' in a continuum 'X' is a 'dispersion point' if the remaining
# space, X \ {x}, is 'totally disconnected'.
# A space is 'totally disconnected' if its only connected subsets are single points.
# The problem is to find the maximum possible size of the set of dispersion points.

# Step 2: Formulate the proof that the number of dispersion points is at most 1.
# We will use proof by contradiction.
#
# Hypothesis: Assume a continuum X has two distinct dispersion points, p and q.

# Proof:
# 1. Since q is a dispersion point, the space X \ {q} is totally disconnected.
#
# 2. The point p is an element of X \ {q}. In a totally disconnected space, the
#    connected component of any point is just the point itself.
#
# 3. A known result from topology states that if the component of a point 'p' in a
#    compact metric space Y is just {p}, then for any neighborhood of p, there
#    is a set V containing p that is both open and closed (clopen) in Y.
#    Applying this to our situation (Y = X \ {q}), we can find an open set K
#    in the original space X such that p is in K, and the boundary of K in X
#    is exactly the set {q}.
#
# 4. A component of an open set in a metric space is itself open. We can choose K
#    to be the connected component containing p. Therefore, K is a connected open set
#    in X containing p, with boundary(K) = {q}.
#
# 5. Now, we use our second assumption: p is also a dispersion point. This means
#    that the space X \ {p} is totally disconnected.
#
# 6. Consider the set K from step 4. It is a connected set. Since its boundary is {q},
#    the point p cannot be in K. Thus, K is a subset of X \ {p}.
#
# 7. This leads to a contradiction. K is a connected set (and since it's an open
#    set in a continuum, it contains more than one point). But K is also a subset
#    of the totally disconnected space X \ {p}. A connected subset of a totally
#    disconnected space can have at most one point.
#
# 8. The contradiction arises from our initial hypothesis that two distinct
#    dispersion points, p and q, exist. Therefore, this hypothesis must be false.
#
# Conclusion of proof: A continuum can have at most one dispersion point.

# Step 3: Determine if the maximum value of 1 is achievable.
# The proof shows the number is <= 1. To show the maximum is exactly 1,
# we need to know if a continuum with one dispersion point can exist.
#
# Such spaces do exist. A well-known (though non-trivial) example is a variation
# of the Knaster-Kuratowski fan. The existence of such a space confirms that a
# cardinality of 1 is possible.

# Step 4: Final Answer.
# The maximum cardinality is at most 1, and it is known to be at least 1.
# Therefore, the maximum cardinality is exactly 1.
# The final equation is: max_cardinality = 1.

final_answer = 1
print(final_answer)