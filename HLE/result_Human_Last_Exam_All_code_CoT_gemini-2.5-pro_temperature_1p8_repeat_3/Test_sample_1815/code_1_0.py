# The problem is a mathematical question about topology, not a computational one.
# Based on the mathematical reasoning, there are no such topologies.
# The number of topologies satisfying both conditions is 0.

# Let's represent the final numerical answer.
final_answer = 0

print(final_answer)

# The argument in steps:
# 1. Let tau be a group topology on Z. Since Z is countable, tau is first-countable.
# 2. The condition "no nontrivial convergent sequences" on a first-countable space
#    implies that there are no limit points, which means the topology must be discrete
#    if we assume the T1 separation axiom (implied by T2/Hausdorff).
# 3. For a general (not necessarily T2) topological group G, the condition implies that
#    the quotient group G/cl({e}) must be discrete. For G=Z, cl({0}) = K, which
#    is a subgroup kZ for some k>=0.
# 4. G/K being discrete means K must be an open set.
# 5. The "no nontrivial convergent sequences" condition further implies that the
#    subspace topology on K must be such that 0 is not a limit point of K\{0}.
# 6. This means there is a neighborhood U of 0 such that U intersect K = {0}.
# 7. Since K itself is also a neighborhood of 0 (from step 4), their intersection
#    U intersect K = {0} must also be a neighborhood of 0.
# 8. This forces the topology to be the discrete topology (where {0} is open).
# 9. The second condition is that the topology must be totally bounded.
# 10. The discrete topology on Z is not totally bounded because the neighborhood {0}
#     requires infinitely many translates to cover Z.
# 11. Since the only topology satisfying the first condition (no nontrivial convergent sequences)
#     does not satisfy the second condition (totally bounded), there are no such topologies.
# 12. Hence, the number is 0.