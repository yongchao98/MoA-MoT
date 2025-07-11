import numpy as np
from scipy.optimize import fsolve

def solve_for_p():
    """
    Solves for the mixed strategy equilibrium probability p for N=8.
    """
    N = 8.0
    M = 24.0

    # Define the function whose root we want to find.
    # f(p) = p - g(p) = 0  is equivalent to p = g(p)
    def f(p):
        if p >= 1.0 or p <= 0:  # p must be a probability between 0 and 1
            # Add a penalty for out-of-bounds values to guide the solver
            return 1e6
        
        X = (1 - p / N) ** M
        Z = (p * (N - 1) / N) ** M
        
        # Avoid division by zero if Z is close to 1
        if 1 - Z < 1e-9:
            return 1e6

        return p - (1 - X) / (1 - Z)

    # Use an initial guess close to 1, as analysis suggests p is high.
    p_initial_guess = 0.95
    p_solution = fsolve(f, p_initial_guess)[0]

    return p_solution

def main():
    """
    Main function to calculate and print the final answer.
    """
    p = solve_for_p()
    
    # Calculate the value required by the problem.
    result = np.floor(10000 * (1 - p))
    
    # Final answer needs to output the equation with numbers plugged in
    # This is a bit unusual, the code would solve it, not print the formula.
    # To meet the output requirement, I will manually construct the print statement
    # with the solved value of p.
    N = 8
    M = 24
    print(f"For N = {N} and M = {M}, the mixed strategy equilibrium involves playing")
    print(f"the discrete strategy with probability p.")
    print("The equation for p is derived from the equality of payoffs E_D = E_S_N:")
    print(f"p = (1 - (1 - p/{N})^ {M}) / (1 - (p*({N}-1)/{N})^{M})")
    
    # We substitute the found p to show the equation holds
    p_val = p
    rhs = (1-(1-p_val/N)**M) / (1-(p_val*(N-1)/N)**M)
    
    print(f"\nSolving this numerically gives p ≈ {p_val:.6f}")
    print(f"Let's check: {p_val:.6f} ≈ (1 - (1 - {p_val:.6f}/{N})^{M}) / (1 - ({p_val:.6f}*({N-1})/{N})^{M})")
    print(f"{p_val:.6f} ≈ (1 - {(1-p_val/N):.6f}^{M}) / (1 - {(p_val*(N-1)/N):.6f}^{M})")
    print(f"{p_val:.6f} ≈ (1 - {(1-p_val/N)**M:.6f}) / (1 - {(p_val*(N-1)/N)**M:.6f})")
    print(f"{p_val:.6f} ≈ {(1-(1-p_val/N)**M):.6f} / {(1-(p_val*(N-1)/N)**M):.6f}")
    print(f"{p_val:.6f} ≈ {rhs:.6f}\n")


    print(f"The value to find is floor(10000 * (1-p)).")
    print(f"1 - p ≈ 1 - {p_val:.6f} = {1-p_val:.6f}")
    print(f"10000 * (1 - p) ≈ 10000 * {1-p_val:.6f} = {10000*(1-p_val):.6f}")
    print(f"floor({10000*(1-p_val):.6f}) = {int(result)}")
    print(f"\nFinal calculated value: {int(result)}")

if __name__ == '__main__':
    main()
