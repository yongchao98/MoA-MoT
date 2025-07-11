import numpy as np
from scipy.optimize import newton

def solve_for_equilibrium_probability():
    """
    This function solves for the mixed strategy equilibrium probability p
    and calculates the final required value.
    """
    N = 8
    
    # The equilibrium condition simplifies to the polynomial equation:
    # x^24 - 8x + 7 = 0, where x = 1 - p/8.
    # We need to solve this for x.
    def equation(x):
        return x**24 - 8 * x + 7

    # We can also provide the derivative to the solver to improve efficiency.
    def equation_prime(x):
        return 24 * x**23 - 8
    
    # We are looking for a root between 0 and 1. An initial guess of 0.9 is reasonable.
    # The other root is x=1 (p=0), which we are not interested in.
    try:
        x_sol = newton(equation, 0.9, fprime=equation_prime)
    except RuntimeError:
        print("Numerical solver failed to converge. Please check the model.")
        return

    # From x, we calculate p.
    # p = 8 * (1 - x)
    p = N * (1 - x_sol)

    # The problem asks for floor(10000 * (1-p))
    result = 10000 * (1 - p)
    
    final_answer = int(np.floor(result))

    # The problem asks for the output from a print function.
    print("The final equation to solve is x^24 - 8*x + 7 = 0, where x = 1 - p/8.")
    print(f"The numerical solution for x is approximately {x_sol:.6f}.")
    print(f"This corresponds to a probability p of approximately {p:.6f}.")
    print(f"The value of (1-p) is approximately {1-p:.6f}.")
    print(f"The value of 10000 * (1-p) is approximately {result:.2f}.")
    print(f"The final answer, floor(10000 * (1-p)), is:")
    print(final_answer)

solve_for_equilibrium_probability()