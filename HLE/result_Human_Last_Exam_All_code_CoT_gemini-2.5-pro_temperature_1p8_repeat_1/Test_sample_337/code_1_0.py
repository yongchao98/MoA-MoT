def solve_hat_puzzle():
    """
    Calculates the maximal probability of release for the 16 prisoners.
    """
    # Number of prisoners in the game
    N = 16

    # The strategy involves a subgroup of n=15 prisoners playing a subgame
    n = 15
    # The subgame strategy is based on a perfect Hamming code, where r=4
    # because n = 2^r - 1
    r = 4

    # Total number of possible hat configurations for N prisoners
    total_configs = 2**N
    print(f"There are 16 prisoners, so there are 2^{N} = {total_configs} total possible hat configurations.")
    
    # In the n=15 subgame, the losing configurations are the "codewords"
    # The number of codewords in a perfect Hamming code is 2^(n-r)
    num_codewords = 2**(n - r)
    print(f"A subgroup of 15 prisoners plays a subgame. The number of 'losing' configurations (codewords) in this subgame is 2^({n}-{r}) = {num_codewords}.")

    # Number of subgame configurations where the first 15 prisoners would win
    # Their strategy wins for any configuration that is NOT a codeword.
    # The hat of the 16th prisoner can be black or white (x2).
    winning_configs_main = (2**n - num_codewords) * 2
    print(f"The number of 16-hat configurations where the subgame succeeds is ({2**n} - {num_codewords}) * 2 = {winning_configs_main}.")

    # Number of configurations where the subgame would fail (all pass).
    # This happens when the first 15 hats form a codeword.
    # The 16th prisoner now has to guess. By guessing a fixed color (e.g., 'White'),
    # they will be correct for exactly half of these cases.
    winning_configs_saved_by_16th = num_codewords * 1
    print(f"The 16th prisoner can 'save' half of the subgame's losing configurations. The number of saved configurations is {num_codewords} * 1 = {winning_configs_saved_by_16th}.")
    
    # The total number of winning configurations is the sum of the two cases
    total_winning_configs = winning_configs_main + winning_configs_saved_by_16th
    
    # The maximal probability is the ratio of winning configurations to total configurations
    probability_numerator = total_winning_configs
    probability_denominator = total_configs
    
    # Simplify the fraction by dividing both by their greatest common divisor, which is num_codewords
    simplified_numerator = probability_numerator // num_codewords
    simplified_denominator = probability_denominator // num_codewords

    print(f"\nThe total number of winning configurations is {winning_configs_main} + {winning_configs_saved_by_16th} = {total_winning_configs}.")
    print(f"The maximal probability is the ratio of winning to total configurations.")
    print(f"Probability = {total_winning_configs} / {total_configs}")
    print(f"This simplifies to the fraction: {simplified_numerator}/{simplified_denominator}")
    
    final_prob = probability_numerator / probability_denominator
    print(f"The decimal value is: {final_prob}")

solve_hat_puzzle()
<<<31/32>>>