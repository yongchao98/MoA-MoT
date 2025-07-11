import math
from scipy.optimize import root_scalar

def solve_for_p_and_calculate():
    """
    Solves the equilibrium equation for p and calculates the final requested value.
    """
    N = 8
    T = 3 * N

    # We need to find the root of the function F(p) = 0, derived from the
    # equilibrium condition Payoff(Discrete) = Payoff(Alternative).
    # F(p) = p - p * (p * (N - 1) / N)**T + (1 - p / N)**T - 1
    def F(p):
        # This function is valid for p in (0, 1].
        if p <= 0 or p > 1:
            return float('inf')
        
        # Note: (p*(N-1)/N) is a single term raised to the power T
        term_B = (p * (N - 1) / N) ** T
        term_A = (1 - p / N) ** T
        
        return p - p * term_B + term_A - 1

    # We use a numerical solver to find the root p. We know from analysis
    # that there are trivial roots at p=0 and p=1, and we seek the non-trivial
    # root which lies between them, close to 1.
    # We search in a bracket like [0.9, 0.999] to find the correct root.
    try:
        sol = root_scalar(F, bracket=[0.9, 0.9999], method='brentq')
        p_star = sol.root
    except (ImportError, ValueError):
        # Fallback to a simple bisection method if scipy is not available.
        low = 0.9
        high = 0.9999
        for _ in range(100):
            mid = (low + high) / 2
            if F(mid) < 0:
                low = mid
            else:
                high = mid
        p_star = (low + high) / 2

    # Calculate the value requested by the user
    value = 10000 * (1 - p_star)
    
    # The final answer is the floor of this value
    result = math.floor(value)

    print(f"For N=8, we solve the equilibrium equation:")
    print(f"p - p * (p * 7/8)^24 + (1 - p/8)^24 - 1 = 0")
    print(f"The equilibrium probability found is p = {p_star:.6f}")
    print(f"The value of 1-p is {1-p_star:.6f}")
    print(f"Then we calculate 10000 * (1-p) = {value:.6f}")
    print(f"The final result, floor(10000 * (1-p)), is: {result}")

solve_for_p_and_calculate()