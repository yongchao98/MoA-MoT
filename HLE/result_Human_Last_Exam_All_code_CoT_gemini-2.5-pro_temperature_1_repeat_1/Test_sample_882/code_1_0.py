import math

def solve_equilibrium_probability():
    """
    This function solves the game theory problem to find the value of 
    floor(10000 * (1-p)).
    """
    
    # Constants from the problem description
    N = 8
    TOTAL_PLAYERS = 3 * N

    def equilibrium_equation(p):
        """
        The equilibrium equation is derived by setting Payoff(Discrete) = Payoff(Uniform).
        After simplification, it becomes f(p) = 0, where f(p) is defined below.
        f(p) = (1 - p/N)**TOTAL_PLAYERS - p * (((N-1)*p)/N)**TOTAL_PLAYERS - (1-p)
        """
        # Trivial roots exist at p=0 and p=1. We seek the internal solution.
        if p <= 0 or p >= 1:
            return 0.0
        
        # Use floating point numbers for calculations
        term1 = (1.0 - p / N)**TOTAL_PLAYERS
        term2 = p * (((N - 1.0) * p) / N)**TOTAL_PLAYERS
        term3 = 1.0 - p
        
        return term1 - term2 - term3

    # Use the bisection method to find the root numerically.
    # From analysis, we know a non-trivial root exists between 0 and 1, specifically
    # close to 1. The interval [0.9, 0.9999] safely brackets the root.
    low = 0.9
    high = 0.9999
    p_solution = 0.0
    
    # Ensure the root is bracketed by the chosen interval
    if equilibrium_equation(low) * equilibrium_equation(high) >= 0:
        print("Error: Could not find a bracketing interval for the root.")
        return

    # Perform 100 iterations for high precision
    for _ in range(100):
        mid = (low + high) / 2.0
        if equilibrium_equation(mid) == 0.0:
            p_solution = mid
            break
        elif equilibrium_equation(mid) * equilibrium_equation(low) < 0:
            high = mid
        else:
            low = mid
    
    p_solution = (low + high) / 2.0

    # As requested, output the numbers in the final equation.
    # The equation is: (1 - p/8)**24 - p * ((7*p)/8)**24 - (1-p) = 0
    print("The numbers defining the equilibrium equation for N=8 are:")
    print(f"Number of races (N): {N}")
    print(f"Total players (3N): {TOTAL_PLAYERS}")
    print(f"Base for discrete-player-avoidance term (N-1): {N-1}")

    # Calculate and print the final result as per the problem's request.
    final_calculation = 10000.0 * (1.0 - p_solution)
    final_answer = math.floor(final_calculation)
    
    print(f"\nThe value of p is approximately: {p_solution:.6f}")
    print(f"The final calculated integer value is: {final_answer}")

solve_equilibrium_probability()
<<<32>>>