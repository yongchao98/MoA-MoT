import numpy as np
from scipy.optimize import fsolve

def solve_equilibrium_p():
    """
    Solves for the mixed strategy equilibrium probability p.
    """
    N = 8
    M = 3 * N

    # The equilibrium equation is derived from the indifference condition
    # u_D(p) = u_D'(p), where D is the discrete strategy and D' is a
    # slight deviation.
    # Let y = p/N. The equation becomes:
    # (1 - (1-y)^M) / (M*y) = (1-y)^(M-1)
    # We define a function f(y) whose root we want to find.
    def f(y):
        # We need to handle the case y=0 to avoid division by zero.
        # At y=0, p=0. LHS is 1 by L'Hopital's rule. RHS is 1.
        # So y=0 is a trivial solution we want to avoid.
        if y == 0:
            return 0
        lhs = (1 - (1 - y)**M) / (M * y)
        rhs = (1 - y)**(M - 1)
        return lhs - rhs

    # We need to solve f(y) = 0 for y.
    # A reasonable initial guess for y=p/N would be around 1/N.
    # Since p should be close to 1, y should be close to 1/N.
    initial_guess_y = 1.0 / N
    # Use a numerical solver to find the root y.
    y_solution = fsolve(f, initial_guess_y)[0]

    # Calculate p from y
    p = N * y_solution

    # Calculate the final value as requested by the user
    # to six significant digits.
    one_minus_p = 1 - p
    result = np.floor(10000 * one_minus_p)
    
    # We print the components to show the calculation.
    # Since we need to output the equation, let's represent it
    # with the solved value of p.
    q_val = 1 - p/N
    print(f"The equilibrium condition is derived from equating payoffs:")
    print(f"Payoff(Discrete) = Payoff(Deviate)")
    print(f"This gives the equation: (1-(1-p/N)^M)/(M*p) = (1-p/N)^(M-1)")
    print(f"For N=8 and M=24, with p ~ {p:.6f}:")
    print(f"(1-(1-{p:.6f}/8)^24)/(24*{p:.6f}) = (1-{p:.6f}/8)^(24-1)")
    
    payoff_d = (1-(1-p/N)**M) / (M*p/N * N)
    payoff_d_prime = (1-p/N)**(M-1)
    print(f"{payoff_d:.6f} = {payoff_d_prime:.6f}")
    
    print("\nCalculated value of p: {:.6f}".format(p))
    print("Calculated value of (1-p): {:.6f}".format(one_minus_p))
    print("Value of 10000 * (1-p): {:.6f}".format(10000 * one_minus_p))
    print("Final answer: floor(10000 * (1-p))")
    print(int(result))

solve_equilibrium_p()