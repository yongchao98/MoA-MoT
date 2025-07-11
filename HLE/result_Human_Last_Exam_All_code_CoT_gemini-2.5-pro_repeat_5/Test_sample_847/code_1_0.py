import numpy as np
from scipy.optimize import minimize_scalar

def solve():
    """
    Calculates the minimal cost per bit of information to sort an array
    based on the given problem's cost structure.
    """

    # Define the cost-per-bit function for a general yes/no question
    # where 'yes' has probability p.
    # Expected Cost C(p) = p*1 + (1-p)*3 = 3 - 2p
    # Expected Information H(p) = -p*log2(p) - (1-p)*log2(1-p)
    # Cost per bit f(p) = C(p) / H(p)
    def cost_per_bit(p):
        # p must be in (0, 1) for entropy to be defined and non-zero.
        if p <= 1e-9 or p >= 1 - 1e-9:
            return float('inf')
        
        entropy = -p * np.log2(p) - (1 - p) * np.log2(1 - p)
        expected_cost = 3 - 2 * p
        
        return expected_cost / entropy

    # Find the minimum of the cost_per_bit function for p in (0, 1)
    # using a numerical optimization method.
    result = minimize_scalar(cost_per_bit, bounds=(0, 1), method='bounded')

    min_cost_per_bit_general = result.fun
    optimal_p = result.x
    cost_per_bit_comparison = 2.0

    # The optimal strategy is the one with the minimum cost per bit.
    final_min_cost = min(min_cost_per_bit_general, cost_per_bit_comparison)

    print("To determine the permutation of n elements, we need to acquire log2(n!) bits of information.")
    print("\nThere are two types of questions to acquire this information:")
    print(f"1. Comparison questions: These provide 1 bit of information for a cost of 2 coins. The cost per bit is {cost_per_bit_comparison:.3f}.")
    print(f"2. General questions: A question with a 'yes' probability 'p' has an expected cost of (p * 1) + ((1-p) * 3) and provides H(p) = -p*log2(p) - (1-p)*log2(1-p) bits of information.")
    
    print("\nBy optimizing the probability 'p', we can find the minimum cost per bit for general questions.")
    print(f"The optimal 'yes' probability is found to be p = {optimal_p:.3f}.")
    print(f"This leads to a minimal cost per bit of {min_cost_per_bit_general:.3f} for general questions.")

    print(f"\nComparing the two options ({min_cost_per_bit_general:.3f} vs {cost_per_bit_comparison:.3f}), the general question strategy is more efficient.")
    
    print("\nThe final equation for the total cost to sort the array for large n is:")
    print(f"Total Cost ≈ {final_min_cost:.3f} * n * log2(n)")
    print("The question asks for the constant factor in this expression.")
    
    # Final answer formatting
    print(f"\n<<< {final_min_cost:.3f} >>>")

solve()