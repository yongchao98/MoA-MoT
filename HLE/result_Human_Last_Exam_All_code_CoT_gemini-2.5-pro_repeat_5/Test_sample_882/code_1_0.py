import math

def solve_game_theory_problem():
    """
    Solves the game theory problem to find the mixed strategy equilibrium probability p
    and then calculates the final requested value.
    """
    N = 8
    M = 3 * N - 1  # Number of other players (23)
    C = 3 * N      # Total number of players in some formulas (24)

    # Payoff function for the Discrete strategy
    def u_discrete(p):
        if p == 0:
            # Using L'Hopital's rule for the limit as p -> 0
            return 1.0
        # The general formula for p > 0
        return (1 - (1 - p / N)**C) / (3 * p)

    # Payoff function for the Split-N strategy
    def u_split(p):
        total_payoff = 0
        for k in range(1, N + 1):
            term = math.comb(N, k) * (1 - p * k / N)**M
            if k % 2 == 1:  # if k is odd
                total_payoff += term
            else:  # if k is even
                total_payoff -= term
        return total_payoff

    # We need to solve u_discrete(p) = u_split(p), which is equivalent to finding the root of h(p) = 0
    def h(p):
        return u_split(p) - u_discrete(p)

    # Bisection method to find the root of h(p)
    # Analysis shows the stable equilibrium is the larger root in (0, 1).
    # We can search in a smaller bracket based on initial exploration.
    low = 0.8
    high = 0.9
    
    # Check if the root is within the bracket
    if h(low) * h(high) >= 0:
        # Fallback to a wider search if our initial guess is wrong
        low, high = 0.0, 1.0
        if h(low) * h(high) >= 0:
             print("Could not find a root in the interval.")
             return

    for _ in range(100):  # 100 iterations for high precision
        mid = (low + high) / 2
        if h(mid) < 0:
            low = mid
        else:
            high = mid
    
    p = (low + high) / 2

    # Calculate the final result
    one_minus_p = 1 - p
    result_multiplied = 10000 * one_minus_p
    final_answer = math.floor(result_multiplied)

    # Output the steps as requested
    print(f"The equilibrium probability p is the stable root of U_Discrete(p) = U_Split(p).")
    print(f"Using a numerical solver for N=8, we find p.")
    print(f"p = {p:.6f}")
    print("\nNow we calculate the final value:")
    print(f"1 - p = {one_minus_p:.6f}")
    print(f"10000 * (1 - p) = {result_multiplied:.6f}")
    print(f"floor(10000 * (1 - p)) = {final_answer}")

solve_game_theory_problem()