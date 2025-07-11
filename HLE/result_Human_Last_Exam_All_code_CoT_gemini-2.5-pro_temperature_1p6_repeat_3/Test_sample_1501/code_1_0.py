import sympy

def solve_sqs_questions():
    """
    This function formulates the answers to the user's questions
    based on the combinatorial analysis of the SQS doubling construction.
    """

    # Part (a): The statement is taken as the definition of an ND-pair.
    ans_a = "True"

    # Part (b): Formulate the expression for the new multiplicity.
    # The new multiplicity is the old multiplicity (mu) plus an
    # additional term derived from block counting, which is v/2.
    v, mu = sympy.symbols('v mu')
    b_expression = mu + v / 2
    
    # We create the string representation for the final output.
    ans_b_str = "mu + v/2"

    # Part (c): Based on average multiplicity calculation, it's not a necessity.
    ans_c = "No"

    # Construct the final answer string in the specified format.
    final_answer_string = f"(a) {ans_a}; (b) {ans_b_str}; (c) {ans_c}"

    print(final_answer_string)

solve_sqs_questions()