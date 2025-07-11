def solve_sqs_problem():
    """
    This function provides the solution to the SQS doubling construction problem.

    The solution is based on the following interpretations:
    1. An 'ND-pair' in the SQS(2v) is a pair of points from the same 'half'
       of the construction, i.e., of the form {(x, i), (y, i)}.
    2. The 'multiplicity' of a pair is its lambda_2 value, i.e., the number
       of blocks that contain it.
    3. The doubling construction is the standard one for SQS(v) -> SQS(2v)
       using a 1-factorization of K_v.
    """

    # Part (a) Solution:
    # A point (x, i) is contained in ND-pairs of the form {(x, i), (y, i)}
    # for all y in Q where y != x. There are v-1 such points y.
    answer_a = "True"

    # Part (b) Solution:
    # Let mu be the multiplicity in SQS(v). So, mu = (v-2)/2.
    # The new multiplicity mu' in SQS(2v) is v-1.
    # We express mu' in terms of mu.
    # From mu = (v-2)/2, we get v = 2*mu + 2.
    # So, mu' = v - 1 = (2*mu + 2) - 1 = 2*mu + 1.
    # The final equation is mu' = 2 * mu + 1.
    # The numbers in this equation are the coefficient 2 and the constant 1.
    coefficient = 2
    constant = 1
    answer_b = f"{coefficient} * mu + {constant}"

    # Part (c) Solution:
    # The multiplicity of any ND-pair in the new SQS(2v) is v-1.
    # For a multiplicity of v, we would need v-1 = v, which implies -1 = 0.
    # This is impossible.
    answer_c = "No"

    # Print the final answer in the required format.
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)


solve_sqs_problem()
<<< (a) True; (b) 2 * mu + 1; (c) No >>>