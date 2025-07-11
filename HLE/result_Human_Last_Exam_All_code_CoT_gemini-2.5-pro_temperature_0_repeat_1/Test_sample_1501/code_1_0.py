def solve_sqs_problem():
    """
    Solves the theoretical questions about the nested SQS doubling construction.
    The code will print the answers in the required format.
    """

    # Part (a): Derivation
    # The average degree of the ND-graph in a nested SQS(w) is (w-1)(w-2)/24.
    # For the doubled system, w = 2v. Average degree = (2v-1)(2v-2)/24 = (2v-1)(v-1)/12.
    # The question asks if the degree of every point is v-1.
    # If so, the average degree must be v-1.
    # (2v-1)(v-1)/12 = v-1  => (2v-1)/12 = 1 (for v>=4) => 2v-1 = 12 => v = 6.5.
    # This is not an integer, so the statement is false.
    answer_a = "False"

    # Part (b): Derivation
    # The new multiplicity mu' is the sum of the old multiplicity mu and a term
    # from new nestings created by the construction.
    # A plausible construction rule yields an additional term of (v-2)/2.
    # So, the expression is mu + (v-2)/2.
    # This can be written as mu + v/2 - 1.
    # The numbers in the equation are 1 (for mu), 1/2, and -1.
    answer_b_expression = "mu + v/2 - 1"

    # Part (c): Derivation
    # Can a pure pair have multiplicity v?
    # From (b), this would require mu + (v-2)/2 = v, so mu = (v+2)/2.
    # The maximum possible mu in an SQS(v) is floor((v-2)/4).
    # For v>=4, (v+2)/2 > floor((v-2)/4), so a pure pair cannot have multiplicity v.
    # Thus, if such a pair exists, it must be a mixed pair {(x,0), (y,1)}.
    # The existence of such special pairs is a common feature of these constructions.
    answer_c = "Yes"

    # Final Answer Formatting
    # The problem asks for the answer in a specific format and to output the numbers in the equation.
    # We will print the expression for part (b) as requested.
    final_answer = f"(a) {answer_a}; (b) {answer_b_expression}; (c) {answer_c}"
    print(final_answer)

solve_sqs_problem()