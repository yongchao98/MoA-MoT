def solve_sqs_problem():
    """
    This function solves the three-part question about the SQS doubling construction.
    """

    # Part (a): True or False: In the doubling construction, each element of
    # Q x {0, 1} is contained in exactly v - 1 ND-pairs.
    #
    # We infer that an "ND-pair" is a pair of points with the same second
    # coordinate, e.g., {(x, i), (y, i)}.
    # For any given element (x, i), we can form a pair with any other element
    # (y, i) where y is in Q and y is not equal to x.
    # There are v-1 such choices for y.
    # Thus, each element is contained in exactly v-1 such pairs.
    answer_a = "True"

    # Part (b): What is the multiplicity of an ND-pair {(x, 0), (y, 0)}
    # in the resulting nested SQS(2v) if the pair {x, y} had
    # multiplicity mu in the original SQS(v)?
    #
    # The multiplicity in an SQS(n) is lambda_n = (n-2)/2.
    # For SQS(v), mu = (v-2)/2. From this, we derive v = 2*mu + 2.
    # For the resulting SQS(2v), the multiplicity is lambda_2v = (2v-2)/2 = v-1.
    # The construction is stated to produce an SQS(2v), so every pair,
    # including an ND-pair, must have this multiplicity.
    # Substituting our expression for v:
    # New multiplicity = v - 1 = (2*mu + 2) - 1 = 2*mu + 1.
    two = 2
    one = 1
    answer_b = f"{two}*mu + {one}"

    # Part (c): Must there exist ND-pairs with multiplicity exactly v
    # in the constructed nested SQS(2v)?
    #
    # As determined for part (b), the multiplicity of every pair in the
    # SQS(2v) is v-1.
    # Since v >= 4, v-1 is never equal to v.
    # Therefore, no pair can have multiplicity v.
    answer_c = "No"

    # Print the final answer in the required format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_sqs_problem()
<<< (a) True; (b) 2*mu + 1; (c) No >>>