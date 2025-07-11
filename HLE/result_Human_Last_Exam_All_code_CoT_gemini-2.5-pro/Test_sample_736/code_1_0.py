import numpy as np
from scipy.optimize import fsolve

def find_maximum_c():
    """
    This function finds the maximum real number c by solving the problem's continuous analogue.

    The constant c is the infimum of the minimum eigenvalues of the matrices A_n.
    This infimum is the limit as n -> infinity, which can be found by solving an
    eigenvalue problem for a continuous integral operator. This leads to a
    transcendental equation z*tanh(z) = 1. The constant c is then given by 1/(2*z^2),
    where z is the solution to the equation.
    """

    # We define the function whose root we need to find.
    # The equation is z * tanh(z) - 1 = 0.
    def equation_to_solve(z):
        return z * np.tanh(z) - 1

    # An initial guess for the root. Based on plotting, the root is near 1.2.
    initial_guess = 1.2

    # Use a numerical solver to find the positive root of the equation.
    z_solution = fsolve(equation_to_solve, initial_guess)[0]

    # Calculate the constant c from the solution z.
    c = 1 / (2 * z_solution**2)

    print("The problem reduces to finding the smallest eigenvalue of a continuous operator.")
    print("This requires solving the transcendental equation: z * tanh(z) = 1")
    print(f"The numerical solution for z is approximately: {z_solution:.8f}")
    print("The maximum value of c is then calculated as c = 1 / (2 * z^2).")
    print(f"The final numerical value for c is approximately: {c:.8f}")

if __name__ == '__main__':
    find_maximum_c()