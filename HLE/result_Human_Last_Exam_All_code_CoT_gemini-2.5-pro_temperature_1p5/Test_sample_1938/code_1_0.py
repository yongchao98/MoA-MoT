import math

def solve():
    """
    Solves for the minimum q such that Alice's win probability is >= 95%.
    """

    # This function calculates P_A for a given q based on the derived formula.
    # P_A(q) = ( 3q^2 - sqrt(4 - 3q^4) ) / (2q^6)
    def calculate_P_A(q):
        # For P_A to be real and non-negative, q^4 must be >= 1/3.
        # Below this critical value, Alice's win probability is 0.
        if q**4 < 1/3:
            return 0.0
        
        # Avoid division by zero, though q will be positive in our search range.
        if q == 0:
            return 0.0

        try:
            numerator = 3 * q**2 - math.sqrt(4 - 3 * q**4)
            denominator = 2 * q**6
            return numerator / denominator
        except ValueError:
            # This can happen if 4 - 3*q^4 is negative, but our q range [0,1] prevents this.
            return 0.0

    # We need to find q_0 such that P_A(q_0) = 0.95.
    # We can use binary search because P_A(q) is monotonically increasing with q.
    
    # The lowest possible q for a non-zero win probability.
    low = (1/3)**0.25 
    # The highest possible q. At q=1, P_A=1.
    high = 1.0
    
    target_prob = 0.95

    # Perform binary search for 100 iterations for high precision.
    for _ in range(100):
        mid_q = (low + high) / 2
        prob_at_mid = calculate_P_A(mid_q)

        if prob_at_mid < target_prob:
            low = mid_q
        else:
            high = mid_q
            
    # The result of the search is our q_0, which is stored in 'high'.
    q0 = high
    
    # The problem asks for floor(100 * q_0).
    result = math.floor(100 * q0)
    
    # We print the final required value.
    # The final equation is result = floor(100 * q0)
    # The numbers in this equation are 100, q0, and result.
    # For example: floor(100 * 0.9549...) = 95
    print(f"The minimum q is approximately: {q0}")
    print(f"The final calculation is: floor(100 * {q0})")
    print(f"The result is: {result}")


solve()