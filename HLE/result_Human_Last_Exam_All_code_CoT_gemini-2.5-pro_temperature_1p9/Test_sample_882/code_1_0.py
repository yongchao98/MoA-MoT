import math

def solve_for_p():
    """
    Solves the equation 1 - p = (1 - p/8)^24 for p using a numerical method.
    The equation describes the mixed strategy Nash Equilibrium condition for N=8.
    """
    from scipy.optimize import fsolve

    N = 8
    
    # Define the function whose root we want to find.
    # f(p) = (1 - p/N)^(3*N) - (1-p)
    func = lambda p: (1 - p/N)**(3*N) - (1 - p)

    # Initial guess for p. Since the discrete strategy was a NE for N<8,
    # p is expected to be close to 1 for N=8.
    p_initial_guess = 0.95
    
    # Use fsolve to find the root.
    p_solution = fsolve(func, p_initial_guess)[0]
    
    return p_solution

def main():
    """
    Main function to execute the plan and print the final result.
    """
    N = 8
    p = solve_for_p()
    one_minus_p = 1 - p
    result = math.floor(10000 * one_minus_p)

    # We are asked to output each number in the final equation.
    # So we demonstrate the final calculation with the solved value of p.
    # The value of p is given to six significant digits.
    p_six_sig_figs = float(f"{p:.6g}")
    one_minus_p_val = 1 - p_six_sig_figs
    final_val = 10000 * one_minus_p_val

    print(f"For N={N}, we must solve the equation: 1 - p = (1 - p/{N})^({3*N})")
    print(f"Solving 1 - p = (1 - p/8)^24 numerically gives p = {p_six_sig_figs}")
    print(f"Then, the value of 1-p is {one_minus_p_val:.6f}")
    print(f"The calculation is: floor(10000 * {one_minus_p_val:.6f}) = floor({final_val:.2f}) = {result}")

if __name__ == '__main__':
    main()
