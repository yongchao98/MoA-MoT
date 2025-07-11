import numpy as np
import math

def calculate_cost_per_bit(p):
    """
    Calculates the expected cost per bit for a Type 2 question.
    
    Args:
        p: The probability of a "yes" answer.

    Returns:
        The cost per bit, or infinity if p is 0 or 1.
    """
    if p <= 1e-9 or p >= 1 - 1e-9:
        return float('inf')
        
    # H(p) is the Shannon entropy for a binary event with probability p
    entropy = -p * math.log2(p) - (1 - p) * math.log2(1 - p)
    
    # C(p) is the expected cost
    expected_cost = 1 * p + 3 * (1 - p)
    
    if entropy == 0:
        return float('inf')
        
    return expected_cost / entropy

def find_minimal_cost():
    """
    Finds the minimal cost per bit by numerically searching for the minimum
    of the cost function.
    """
    # Create a range of probabilities to test, from 0 to 1
    probabilities = np.linspace(0.0001, 0.9999, 200000)
    
    # Calculate the cost per bit for each probability
    costs = np.array([calculate_cost_per_bit(p) for p in probabilities])
    
    # Find the minimum cost
    min_cost = np.min(costs)
    
    return min_cost

if __name__ == '__main__':
    min_cost_per_bit = find_minimal_cost()
    
    # The total cost for sorting an array of size n for large n is given by:
    # Cost(n) ≈ C * n * log2(n)
    # where C is the minimal cost per bit we have found.
    # The problem asks for this coefficient C, rounded to 3 decimal places.
    
    # We output the values that form the final cost equation.
    print(f"The minimal cost is approximately: C * n * log2(n)")
    print(f"The coefficient C, representing the minimal cost per bit, is: {min_cost_per_bit:.3f}")
    # Final answer tag expects only the number.
    # print(f"<<<{min_cost_per_bit:.3f}>>>")
