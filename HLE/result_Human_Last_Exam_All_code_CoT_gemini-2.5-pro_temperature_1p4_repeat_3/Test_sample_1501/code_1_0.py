def solve_sqs_doubling():
    """
    Solves the theoretical questions about the SQS doubling construction.
    """

    # Part (a) analysis
    # Let d(p, S) be the number of ND-pairs in a system S containing point p.
    # In SQS(v), the number of blocks through a point x is r = (v-1)(v-2)/6.
    # Each such block contributes one ND-pair containing x. So, d(x, S) = r.
    # In SQS(2v), for a point (x,0), the ND-pairs containing it are:
    # 1. From Type 1 blocks: d(x, S) = (v-1)(v-2)/6 pairs of the form {(x,0), (y,0)}.
    # 2. From Type 2 blocks: (v-1) pairs of the form {(x,0), (x,1)}.
    # Total count d((x,0), S') = (v-1)(v-2)/6 + (v-1).
    # This simplifies to the equation: d((x,0), S') = (v-1)(v+4)/6.
    # The question is whether d((x,0), S') == v-1.
    # (v-1)(v+4)/6 == v-1  => (v+4)/6 == 1 => v+4 == 6 => v == 2.
    # This is not true for v >= 4.
    answer_a = "False"

    # Part (b) analysis
    # Let m(P, S) be the multiplicity of pair P in system S.
    # We want m({(x,0),(y,0)}, S').
    # These pairs only come from Type 1 nestings. A pair {(x,0),(y,0)} is formed
    # if and only if {x,y} was an ND-pair in S.
    # The construction generates one {(x,0),(y,0)} for each {x,y} in the original nesting.
    # So, the final equation for the new multiplicity m_new is: m_new = mu.
    answer_b_expression = "mu"
    
    # Part (c) analysis
    # We check if any ND-pair can have multiplicity equal to v.
    # - Horizontal pairs {(x,0),(y,0)} have multiplicity mu. mu <= (v-2)/2 < v for v>=4.
    # - Vertical pairs {(x,0),(x,1)} have multiplicity v-1. This is not v.
    # - Other pairs have multiplicity 0.
    # No pair type has multiplicity v.
    answer_c = "No"

    final_answer = f"(a) {answer_a}; (b) {answer_b_expression}; (c) {answer_c}."
    print(final_answer)

solve_sqs_doubling()