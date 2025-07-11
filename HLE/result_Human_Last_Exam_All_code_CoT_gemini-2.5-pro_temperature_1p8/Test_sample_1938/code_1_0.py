import math

def solve_for_q0():
    """
    This function solves for the minimum value of q, q_0, such that Alice's
    winning probability is at least 95%.
    """

    # We derived the following equation for q, by setting Alice's winning
    # probability p_A to 0.95:
    # q * (1 - (1 - p_A) * q)^3 = 1 - (1 - p_A)^(1/3)
    # Substituting p_A = 0.95, the specific equation to solve is:
    # q * (1 - 0.05 * q)^3 = 1 - 0.05^(1/3)

    p_a = 0.95
    c1 = 1.0 - p_a
    c2 = c1**(1.0 / 3.0)
    target_value = 1.0 - c2

    print("The problem reduces to solving the following equation for q:")
    # The user request asks to "output each number in the final equation".
    # The numbers in the equation are: 1, 3, 0.05, 1, 1/3, 0.05
    print(f"q * (1 - {c1:.2f}*q)^3 = 1 - ({c1:.2f})^(1/3)")
    print(f"Numerically, we are solving: q * (1 - {c1:.2f}*q)^3 = {target_value}")
    print("-" * 20)

    # We use binary search to find the value of q in the range [0, 1].
    # The function f(q) = q * (1 - 0.05*q)^3 is monotonically increasing for q in [0,1].
    low_q = 0.0
    high_q = 1.0
    
    # Iterate for 100 times to get a very precise value for q.
    for _ in range(100):
        mid_q = (low_q + high_q) / 2
        
        # Calculate the value of the left side of the equation.
        current_value = mid_q * (1 - c1 * mid_q)**3
        
        if current_value < target_value:
            # If the current value is too small, q needs to be larger.
            low_q = mid_q
        else:
            # If the current value is too large, q can be smaller.
            high_q = mid_q

    # After the search, high_q (or low_q) is a very good approximation of q_0.
    q0 = high_q
    
    # The question asks for the floor of 100 * q_0.
    result = math.floor(100 * q0)

    print(f"The minimum value of q is q_0 ≈ {q0:.6f}")
    print(f"100 * q_0 ≈ {100 * q0:.6f}")
    print("The final integer answer is:")
    print(result)

solve_for_q0()