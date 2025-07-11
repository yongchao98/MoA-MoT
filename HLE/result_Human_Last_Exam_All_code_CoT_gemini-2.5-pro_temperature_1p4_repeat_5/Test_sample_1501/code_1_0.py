def solve_sqs_problem():
    """
    This function solves the questions based on the properties of Steiner Quadruple Systems.

    Plan:
    1. A system that is an SQS(k) has the property that any pair of points is
       contained in exactly lambda = (k-2)/2 blocks. For the resulting SQS(2v),
       the multiplicity of any pair is (2v-2)/2 = v-1.

    2. (a) True or False: Each element of Q x {0, 1} is contained in exactly v - 1 ND-pairs.
       If we define an "ND-pair" as a pair of points from the same partition,
       e.g., {(x, 0), (y, 0)}, then an element (a, 0) is contained in v-1 such pairs:
       {(a, 0), (y, 0)} for all y != a. This makes the statement True.

    3. (b) What is the multiplicity of an ND-pair {(x, 0), (y, 0)}?
       Based on step 1, its multiplicity is v-1.
       We are asked to express this in terms of mu, the multiplicity in the original SQS(v).
       In an SQS(v), mu = (v-2)/2.
       Solving for v: v = 2*mu + 2.
       Substituting into v-1 gives: (2*mu + 2) - 1 = 2*mu + 1.
       The numbers in this equation are 2 and 1.

    4. (c) Must there exist ND-pairs with multiplicity exactly v?
       No. All pairs in an SQS(2v) have multiplicity v-1. Since v >= 4, v != v-1.
    """

    # Part (a)
    answer_a = "True"

    # Part (b)
    # The new multiplicity is 2*mu + 1.
    # The numbers in the equation are 2 and 1.
    mu_symbol = "mu"
    expression_b = f"2*{mu_symbol} + 1"

    # Part (c)
    answer_c = "No"

    # Construct the final output string
    final_answer = f"(a) [{answer_a}]; (b) [{expression_b}]; (c) [{answer_c}]."

    print(final_answer)

solve_sqs_problem()