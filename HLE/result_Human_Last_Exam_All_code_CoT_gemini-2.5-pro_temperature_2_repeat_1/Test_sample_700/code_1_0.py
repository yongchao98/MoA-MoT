# This problem is theoretical and does not require coding.
# The reasoning for the complexity classes is as follows:
# G = {M | exists w forall s: M accepts sw}
# This is equivalent to: exists w, for all reachable states q, delta(q, w) = q_f.
# The structure is exists-forall. A witness for the forall-part's falsity
# (a reachable state q that doesn't go to q_f) can be small (path length < n).
# So the inner forall-part is in co-NP.
# The outer exists-part makes the whole problem reside in Sigma_2^P.
#
# G' = {(k,M) | exists w forall s of length k: M accepts sw}
# This is equivalent to: exists w, for all states q reachable in k steps, delta(q, w) = q_f.
# The set of states reachable in k steps, S_k, can be computed in polynomial time
# in |M| and log k (using matrix exponentiation).
# After computing S_k, checking the condition for a given w is a simple P-time loop.
# The structure is: exists w such that a P-time predicate holds. This is NP.

print("Sigma_2, NP")
