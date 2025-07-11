def solve_quantum_dna_probability():
    """
    Calculates the limiting probability P(n) as n -> infinity for the
    quantum DNA polymerase problem.
    """
    # There are 8 distinct nucleotide bases.
    num_bases = 8

    # For any single base 'x', the collapse condition is met if its count, c_x,
    # satisfies c_x mod 8 = x. In the limit of large n, the distribution of
    # c_x mod 8 is uniform over {0, 1, ..., 7}.
    # The probability of collapse for one base is therefore 1/8.
    
    # The probability of success (no collapse) for a single base is 1 minus
    # the probability of failure.
    # P(success for one base) = 1 - 1/8 = 7/8
    prob_success_one_base_num = num_bases - 1
    prob_success_one_base_den = num_bases

    # For the entire sequence to be replicated, the success condition must hold
    # for all 8 bases simultaneously. In the limit, these events are independent.
    # The total probability is the product of individual probabilities.
    # P(total success) = (P(success for one base)) ^ num_bases = (7/8)^8

    # Calculate the numerator and denominator of the final probability fraction.
    numerator = prob_success_one_base_num ** num_bases
    denominator = prob_success_one_base_den ** num_bases

    # Print the breakdown of the calculation as requested.
    print(f"The closed-form expression for the limiting probability is ({prob_success_one_base_num}/{prob_success_one_base_den})^{num_bases}.")
    print(f"The final equation is derived from (P_success_one_base) ^ (num_bases):")
    print(f"  = ( ( {num_bases} - 1 ) / {num_bases} ) ^ {num_bases}")
    print(f"  = ( {prob_success_one_base_num}^{num_bases} ) / ( {prob_success_one_base_den}^{num_bases} )")
    print(f"\nCalculating each part:")
    print(f"Numerator: {prob_success_one_base_num}^{num_bases} = {numerator}")
    print(f"Denominator: {prob_success_one_base_den}^{num_bases} = {denominator}")
    print(f"\nThe final probability as a fraction is: {numerator} / {denominator}")
    
    # The result is exact because the denominator is a power of 2.
    final_probability_decimal = numerator / denominator
    print(f"As a decimal, the limiting probability is: {final_probability_decimal}")

solve_quantum_dna_probability()