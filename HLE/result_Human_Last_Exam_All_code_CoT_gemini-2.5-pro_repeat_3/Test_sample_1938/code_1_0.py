import math

def solve():
    """
    This function calculates the minimum q such that Alice's winning probability is at least 95%,
    and then finds floor(100*q).
    """

    # This is the function for Alice's winning probability p_A(q).
    # We want to find q such that p_A(q) = 0.95.
    def p_A(q):
        # The non-trivial solution for p_A exists only for q^4 > 1/3
        if q**4 <= 1/3:
            return 0
        
        q2 = q * q
        q4 = q2 * q2
        q6 = q4 * q2
        
        # The equation for Alice's win probability is p_A = (3q^2 - sqrt(4 - 3q^4)) / (2q^6)
        numerator = 3 * q2 - math.sqrt(4 - 3 * q4)
        denominator = 2 * q6
        
        return numerator / denominator

    # We need to solve p_A(q) - 0.95 = 0. We will use bisection search.
    target = 0.95
    
    # Set up the search range for q.
    # The minimum possible q for a non-zero win probability is (1/3)^(1/4) ~ 0.759
    # The maximum is 1.
    low = (1/3)**0.25
    high = 1.0
    
    # We can narrow the search range based on initial checks.
    # p_A(0.90) is just under 0.95 and p_A(0.91) is just over 0.95.
    low = 0.90
    high = 0.91

    # Perform bisection search for q_0
    for _ in range(100):  # 100 iterations for high precision
        mid_q = (low + high) / 2
        if p_A(mid_q) < target:
            low = mid_q
        else:
            high = mid_q
            
    q0 = high  # The result of the search

    # The final step is to calculate floor(100 * q0)
    result = math.floor(100 * q0)
    
    # As requested, outputting the numbers in the final calculation
    print(f"The minimum value of q is approximately q0 = {q0}")
    print(f"The final equation is floor(100 * q0)")
    print(f"Plugging in the value: floor(100 * {q0}) = floor({100 * q0})")
    print(f"The result is: {result}")

solve()