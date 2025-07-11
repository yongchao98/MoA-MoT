def check_potts_fkg_condition(a, b, c, d):
    """
    Checks the Holley criterion for a single edge in the Potts model.
    The positive correlation property holds if this function returns True for all a,b,c,d.

    Args:
        a: spin of configuration xi at vertex u
        b: spin of configuration eta at vertex u
        c: spin of configuration xi at vertex v
        d: spin of configuration eta at vertex v
    """
    # Left-hand side of the inequality
    lhs = (1 if max(a, b) == max(c, d) else 0) + \
          (1 if min(a, b) == min(c, d) else 0)

    # Right-hand side of the inequality
    rhs = (1 if a == c else 0) + \
          (1 if b == d else 0)

    # The condition is LHS >= RHS
    is_satisfied = lhs >= rhs

    print(f"Checking for a={a}, b={b}, c={c}, d={d}:")
    print(f"  Configuration xi: spin u={a}, spin v={c}")
    print(f"  Configuration eta: spin u={b}, spin v={d}")
    print("-" * 20)
    print(f"I(max({a},{b}) = max({c},{d})) + I(min({a},{b}) = min({c},{d})) >= I({a}={c}) + I({b}={d})")
    print(f"    {lhs} >= {rhs}")
    print(f"Condition satisfied: {is_satisfied}\n")
    return is_satisfied

# For q>=3, we can choose distinct values 1, 2, 3.
# This specific choice serves as a counterexample.
a_val = 2
b_val = 1
c_val = 2
d_val = 3

check_potts_fkg_condition(a_val, b_val, c_val, d_val)

# This counterexample shows the property fails for q>=3 on any graph with an edge.
# The question requires the property to hold for all q>=2.
# This implies the graph cannot have any edges.
# For a connected graph, this means it must be a single vertex, so its max degree is 0.
print("The analysis shows that the condition fails for q>=3 if the graph has any edges.")
print("This implies that for the statement to hold for ALL q>=2, the graph must have no edges.")
print("A connected graph with no edges must be a single vertex, which has a maximum degree of 0.")
print("Therefore, the largest possible value for d is 0.")
