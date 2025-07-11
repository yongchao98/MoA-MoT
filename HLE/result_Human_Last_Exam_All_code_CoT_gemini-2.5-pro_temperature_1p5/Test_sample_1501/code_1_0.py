def solve_sqs_problem():
    """
    This function provides the solution to the nested SQS doubling construction problem.
    """
    # Part (a)
    answer_a = "False"
    
    # Part (b)
    # The multiplicity mu in the original SQS(v) relates to the new multiplicity.
    # The new multiplicity is mu + (v-2)/2.
    answer_b_expr = "mu + (v-2)/2"
    
    # Part (c)
    # The multiplicity of any ND-pair P, M_ND(P), must be less than or equal to
    # the total number of blocks containing P, which is lambda_2' = v-1.
    # Since v >= 4, v > v-1. So the multiplicity cannot be v.
    answer_c = "No"
    
    # Print the formatted answer
    print(f"(a) {answer_a}; (b) {answer_b_expr}; (c) {answer_c}.")

solve_sqs_problem()
# The required format is <<<answer content>>>. The python code is a tool to derive it.
# The final answer is constructed from the logic above.
# (a) False
# (b) \mu + (v-2)/2  (I'll represent mu as mu)
# (c) No
final_answer = "(a) False; (b) mu + (v-2)/2; (c) No."