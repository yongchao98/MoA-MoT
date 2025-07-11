import numpy as np
from scipy.optimize import fsolve
import warnings

def solve_alpha():
    """
    This function sets up and attempts to solve the equation for alpha based on the problem description.
    """
    # Step 1: Define the constants from the problem statement.
    T = np.log(10)
    B = 0.5 * 10**20 / (99**2)

    # Step 2: Define the coefficients of the equation for alpha.
    # The derived equation is of the form: K1 * alpha**5 - K2 * alpha**8 = B
    
    # Calculate intermediate constants k2 and k3
    k2 = 1 - np.exp(-2 * T)  # This is 1 - 1/100 = 0.99
    k3 = 1 - np.exp(-3 * T)  # This is 1 - 1/1000 = 0.999

    # From the derivation, the equation C * A**4 / 4 - A**8 / 8 = B is obtained,
    # where C = (3*alpha)/k3 and A = (2*alpha)/k2.
    # Substituting A and C leads to:
    # ( (3*alpha)/k3 * ((2*alpha)/k2)**4 / 4 ) - ( ((2*alpha)/k2)**8 / 8 ) = B
    # ( 12 / (k3 * k2**4) ) * alpha**5 - ( 32 / k2**8 ) * alpha**8 = B
    K1 = 12 / (k3 * k2**4)
    K2 = 32 / (k2**8)

    # Step 3: Print the final equation with its numerical coefficients.
    print("The derived equation for alpha is:")
    print(f"{K1:.4f} * alpha^5 - {K2:.4f} * alpha^8 = {B:.4e}")
    print("-" * 50)

    # Step 4: Analyze the equation before attempting to solve it.
    # Let g(alpha) = K1 * alpha^5 - K2 * alpha^8.
    # We find the maximum value of g(alpha) for alpha > 0.
    # The derivative g'(alpha) = 5 * K1 * alpha^4 - 8 * K2 * alpha^7 = 0
    # gives the value of alpha for which g(alpha) is maximum.
    alpha_for_max_g = (5 * K1 / (8 * K2))**(1/3)
    max_g_value = K1 * alpha_for_max_g**5 - K2 * alpha_for_max_g**8

    print("Analysis of the equation:")
    print(f"The left-hand side (LHS) of the equation is a function that has a maximum value.")
    print(f"The maximum value of the LHS occurs at alpha = {alpha_for_max_g:.4f}")
    print(f"The maximum value of the LHS is approximately {max_g_value:.4f}")
    print(f"The right-hand side (RHS) of the equation is B = {B:.4e}")
    print("\nConclusion: Since the maximum value of the LHS is a small positive number")
    print("while the RHS is an extremely large positive number, there is no real solution for alpha.")
    print("This suggests a likely typo in the problem statement's constants.")
    print("-" * 50)

    # Step 5: Attempt to solve the equation numerically to confirm the analysis.
    def equation_to_solve(alpha):
        # fsolve finds roots of a function, so we look for where f(alpha) = 0.
        return K1 * alpha**5 - K2 * alpha**8 - B

    print("Attempting to find a numerical solution with fsolve...")
    # Suppress the warning fsolve issues when it fails to converge.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # We provide an initial guess. The value that maximizes the LHS is a reasonable choice.
        initial_guess = alpha_for_max_g
        solution, infodict, ier, mesg = fsolve(equation_to_solve, initial_guess, full_output=True)

    # Check if the solver was successful.
    if ier == 1:
        # This part is not expected to be reached.
        print(f"A numerical solution was found: alpha = {solution[0]}")
    else:
        print("The numerical solver failed to find a solution, as expected from the analysis.")
        print(f"Solver message: {mesg}")

if __name__ == '__main__':
    solve_alpha()