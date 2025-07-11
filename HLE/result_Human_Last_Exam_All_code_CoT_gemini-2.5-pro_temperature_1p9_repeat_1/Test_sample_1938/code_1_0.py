import math

def solve_game_theory_problem():
    """
    This function encapsulates the entire logic for solving the problem.
    It calculates the minimum probability q and provides the final answer as requested.
    """
    
    # Step 1: Set up the problem parameters.
    # We are looking for the minimum probability q such that Alice's winning
    # probability, x, is at least 0.95. We will find q for x = 0.95.
    x0 = 0.95

    # Step 2: Formulate the equation.
    # Let x be the probability that Alice wins, assuming it is her turn.
    # Let y be the probability that Alice wins, assuming it is Bob's turn.
    # The relationships are: x = 1 - (1 - q*y)^3 and y = (q*x)^3.
    # Substituting y gives: x = 1 - (1 - q^4 * x^3)^3.
    #
    # We rearrange this to find the non-zero solution for x. This can be expressed
    # as a polynomial. Let Q = q^4. The equation becomes a cubic polynomial in Q for a fixed x:
    # (x^8) * Q^3 - 3*(x^5) * Q^2 + 3*(x^2) * Q - 1 = 0

    # Step 3: Solve the equation for Q.
    # We use numerical methods to solve this cubic equation for Q, with x = 0.95.
    # The coefficients for the polynomial a*Q^3 + b*Q^2 + c*Q + d = 0 are calculated first.
    
    a = x0**8
    b = -3 * x0**5
    c = 3 * x0**2
    d = -1

    def h(Q_val):
        return a*Q_val**3 + b*Q_val**2 + c*Q_val + d

    # We use a binary search algorithm because we can establish that the
    # single real root lies between 0 and 1.
    low = 0.0
    high = 1.0
    # 100 iterations of binary search are sufficient for high precision.
    for _ in range(100):
        mid = (low + high) / 2
        if h(mid) < 0:
            low = mid
        else:
            high = mid
    
    Q0 = (low + high) / 2
    
    # Step 4: Calculate q0 from Q0.
    q0 = Q0**(1/4)
    
    # Step 5: Calculate the final answer.
    answer = math.floor(100 * q0)

    # Step 6: Print the explanation and results as requested.
    print("The relationship between Alice's win probability (x) and the edge probability (q) is given by:")
    print("x = 1 - (1 - q^4 * x^3)^3")
    print("\nWe want to find the minimum q for which x >= 0.95.")
    print("We solve for q when x = 0.95.")
    print("Let Q = q^4. The equation can be rewritten as a cubic polynomial in Q:")
    print("x^8 * Q^3 - 3*x^5 * Q^2 + 3*x^2 * Q - 1 = 0")
    print("\nSubstituting x = 0.95, the numbers in the final equation are:")
    print(f"a (coefficient of Q^3) = x^8 = {a:.6f}")
    print(f"b (coefficient of Q^2) = -3*x^5 = {b:.6f}")
    print(f"c (coefficient of Q)   = 3*x^2 = {c:.6f}")
    print(f"d (constant term)      = -1.000000")
    print("\nThe equation is:")
    print(f"{a:.6f} * Q^3 + {b:.6f} * Q^2 + {c:.6f} * Q - 1.000000 = 0")
    
    print(f"\nSolving this numerically for Q, we find the root Q_0 ≈ {Q0:.6f}")
    print(f"Since Q = q^4, q is the 4th root of Q.")
    print(f"q_0 = (Q_0)^(1/4) ≈ {q0:.6f}")
    
    print("\nThe problem asks for the value of floor(100 * q_0).")
    print(f"100 * q_0 ≈ {100 * q0:.6f}")
    print(f"floor(100 * q_0) = {answer}")

# Execute the solution function.
solve_game_theory_problem()