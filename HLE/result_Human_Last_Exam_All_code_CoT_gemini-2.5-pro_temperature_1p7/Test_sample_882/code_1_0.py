import math

def solve_equilibrium_probability():
    """
    This function solves for the mixed strategy equilibrium probability 'p'
    and calculates the final requested value.
    """
    N = 8
    num_races = N
    num_players = 3 * N

    # We need to solve the equation f(p) = 0 for the non-trivial root, where:
    # f(p) = 1 - (1 - p / num_races)**num_players - p
    def f(p):
        return 1.0 - (1.0 - p / num_races)**num_players - p

    # We use the bisection method to find the root p in (0, 1).
    # We start the search interval slightly above 0 to avoid the trivial root f(0)=0.
    low = 1e-12
    high = 1.0
    
    # A fixed number of iterations is sufficient for high precision.
    for _ in range(100):
        mid = (low + high) / 2
        # If f(mid) is positive, the root is in the upper half.
        if f(mid) > 0:
            low = mid
        # If f(mid) is negative, the root is in the lower half.
        else:
            high = mid
            
    p = (low + high) / 2
    
    # Calculate the final value as requested by the problem.
    final_value = math.floor(10000 * (1 - p))

    # Print the context and the result.
    print("In the mixed strategy Nash equilibrium, the probability 'p' solves the equation:")
    print(f"p = 1 - (1 - p/{num_races})^{num_players}")
    print(f"For N={num_races}, the equation is: p = 1 - (1 - p/{num_races})^{num_players}")
    print(f"The numerically found value for p is approximately: {p:.6f}")
    
    value_to_floor = 10000 * (1 - p)
    print("\nThe final expression to calculate is: floor(10000 * (1 - p))")
    print(f"1 - p ≈ {1-p:.6f}")
    print(f"10000 * (1 - p) ≈ {value_to_floor:.4f}")
    print(f"The floor of the result is: {final_value}")

solve_equilibrium_probability()
<<<427>>>