import numpy as np
from scipy.optimize import root_scalar

def find_optimal_cost_per_bit():
    """
    Calculates the minimal cost per bit of information by solving
    the cubic equation derived from the cost optimization problem.
    """

    # The optimal strategy leads to the cubic equation x^3 + x - 1 = 0,
    # where x is related to the cost per bit 'w' by x = 2^(-1/w).
    def equation(x):
        return x**3 + x - 1

    # Solve for the real root x0 in the interval [0, 1].
    try:
        sol = root_scalar(equation, bracket=[0, 1], method='brentq')
        x0 = sol.root
    except (ImportError, ValueError) as e:
        print("SciPy is required for numerical solving. Using a pre-calculated value.")
        print(f"Error: {e}")
        # Cardano's formula solution for x^3+x-1=0 is ((sqrt(117)+9)/2)^(1/3) - (2/(sqrt(117)+9))^(1/3)
        # For simplicity, we can use a hardcoded value if scipy fails.
        x0 = 0.6823278038280193

    # The equation to solve is x^3 + x - 1 = 0.
    # The coefficients for x^3, x^2, x, and the constant term are:
    c3 = 1
    c2 = 0
    c1 = 1
    c0 = -1
    
    # Per instructions, outputting the numbers in the final equation.
    print("The cost optimization leads to solving the equation:")
    print(f"{c3}*x^3 + {c2}*x^2 + {c1}*x + {c0} = 0")
    print("-" * 30)

    # Calculate the cost per bit 'w' from the root x0.
    # w = -1 / log2(x0) = -ln(2) / ln(x0)
    cost_per_bit = -np.log(2) / np.log(x0)

    print("The minimal number of coins per bit of information is:")
    print(f"{cost_per_bit:.3f}")
    
    return cost_per_bit

if __name__ == '__main__':
    result = find_optimal_cost_per_bit()
    print(f"\n<<<{result:.3f}>>>")
