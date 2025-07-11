import numpy as np
from scipy.optimize import minimize

def solve():
    """
    This script calculates the minimal number of coins required to sort an array of n elements,
    expressed as a constant factor C* in the total asymptotic cost formula C* * n*log2(n).
    """

    # --- Define Costs ---
    # Cost for a general "yes" answer
    cost_yes = 1
    # Cost for a general "no" answer
    cost_no = 3
    
    # --- Function Definitions ---
    
    def shannon_entropy(p):
        """Calculates the Shannon entropy H(p) for a binary event."""
        # Ensure p is within the valid range (0, 1) to avoid math errors.
        if p <= 0 or p >= 1:
            return 0
        return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

    def cost_per_bit(p):
        """
        Calculates the expected cost per bit of information for a general question
        with a "yes" probability of p.
        """
        # The expected cost equation is: p * cost_yes + (1-p) * cost_no
        expected_cost = p * cost_yes + (1 - p) * cost_no
        
        # The information gain is the entropy
        info_gain = shannon_entropy(p)
        
        # To avoid division by zero when p is 0 or 1
        if info_gain == 0:
            return np.inf
            
        return expected_cost / info_gain

    # --- Optimization ---
    
    # We find the minimum cost per bit by optimizing over p.
    # Initial guess for p is 0.7, as a higher probability for the cheaper "yes" answer is beneficial.
    # Bounds for p are set slightly inside (0, 1) for numerical stability.
    result = minimize(cost_per_bit, x0=0.7, bounds=[(1e-9, 1 - 1e-9)])
    
    # The minimal cost per bit is the minimum value found by the optimizer.
    optimal_cost_per_bit = result.fun

    # --- Outputting the Result ---
    
    # The final equation for total cost is approximately C * n*log2(n).
    # The constant C is the minimum of the cost_per_bit function.
    # The function minimized is (p * 1 + (1-p) * 3) / H(p).
    
    print("The final asymptotic cost to sort the array is given by the equation: Cost ≈ C * n * log2(n).")
    print(f"The constant C is the minimal cost per bit of information, found by minimizing the function (p*{cost_yes} + (1-p)*{cost_no}) / H(p).")
    print(f"The minimal cost per bit, C, is approximately: {optimal_cost_per_bit:.3f}")

solve()