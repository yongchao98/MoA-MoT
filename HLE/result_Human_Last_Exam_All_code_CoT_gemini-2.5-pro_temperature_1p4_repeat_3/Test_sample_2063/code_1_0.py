def solve_limiting_probability():
    """
    Calculates the limiting probability P(n) as n -> infinity.

    The problem states that replication fails if for any base x in {0..7},
    the count of that base (c_x) modulo 8 is equal to x.
    P(success) = P(c_0 mod 8 != 0 AND c_1 mod 8 != 1 AND ... AND c_7 mod 8 != 7)

    As n -> infinity, two key things happen:
    1. For any given base x, the distribution of its count modulo 8 becomes uniform.
       P(c_x mod 8 = k) -> 1/8 for any k in {0..7}.
       Therefore, P(c_x mod 8 != x) -> 1 - 1/8 = 7/8.

    2. The events for different bases become independent.
       P(A and B) -> P(A) * P(B).

    Combining these, the total probability is the product of the individual probabilities.
    """

    # The number of distinct bases, which is also the modulus for the counts.
    num_bases = 8

    # For a single base x, the probability of the failure condition (c_x mod 8 = x)
    # approaches 1/num_bases.
    prob_failure_per_base = 1 / num_bases

    # The probability of the success condition (c_x mod 8 != x) is therefore:
    prob_success_per_base = 1 - prob_failure_per_base

    # The number of independent conditions that must all be met is num_bases.
    exponent = num_bases

    # The final limiting probability is the product of the success probabilities for each base.
    limiting_probability = prob_success_per_base ** exponent

    # The components of the final equation are:
    base_numerator = 7
    base_denominator = 8
    
    print(f"The limiting probability is given by the expression: ({base_numerator}/{base_denominator})^{exponent}")
    print(f"The calculated value is: {limiting_probability}")

solve_limiting_probability()