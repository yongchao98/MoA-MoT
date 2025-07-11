import math

def solve_for_q0():
    """
    This function solves for the minimum probability q (q0) such that Alice's
    winning probability is at least 95%.
    """
    # The condition p_A >= 0.95 leads to the following equation for the minimum q, q0,
    # where p_A = 0.95:
    # q0 * (1 - 0.05 * q0)**3 = 1 - (0.05)**(1/3)
    # We solve this numerically.

    # This is the target value for the function g(q) below.
    target = 1 - (0.05)**(1/3)

    # This is the left-hand side of the equation we need to solve.
    def g(q):
        return q * (1 - 0.05 * q)**3

    # We use binary search to find q0 in the interval [0, 1].
    # The function g(q) is monotonically increasing in this interval.
    low = 0.0
    high = 1.0

    # 100 iterations are sufficient for high precision.
    for _ in range(100):
        mid = (low + high) / 2
        if g(mid) < target:
            low = mid
        else:
            high = mid

    q0 = (low + high) / 2
    
    # The problem asks for floor(100 * q0).
    result = math.floor(100 * q0)
    
    # The following print statements show the components of the final calculation as requested.
    print(f"The equation to find the minimum q, denoted q0, is: q0 * (1 - 0.05*q0)^3 = 1 - 0.05^(1/3)")
    print(f"Solving this numerically gives q0 ≈ {q0:.6f}")
    print(f"The final expression to compute is: floor(100 * q0)")
    print(f"Substituting the numbers: floor(100 * {q0:.6f}) = floor({100 * q0:.6f}) = {result}")

solve_for_q0()