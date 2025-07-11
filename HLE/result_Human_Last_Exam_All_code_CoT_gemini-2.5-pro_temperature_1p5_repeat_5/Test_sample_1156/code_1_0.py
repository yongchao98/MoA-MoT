import numpy as np
from scipy.interpolate import interp1d
from scipy.special import zeta
from scipy.optimize import curve_fit

def perron_frobenius_operator(rho_func, y, k_truncate=100):
    """
    Applies the Perron-Frobenius operator L to a density function rho_func.
    L[rho](y) = 2 * sum_{k=1 to inf} rho((y+k)^-2) / (y+k)^3
    """
    s = 0.0
    for k in range(1, k_truncate + 1):
        x_val = 1.0 / ((y + k) ** 2)
        s += rho_func(x_val) / ((y + k) ** 3)
    return 2.0 * s

def solve_invariant_density(n_points=200, n_iter=10):
    """
    Numerically finds the invariant density by iterating the P-F operator.
    """
    # Grid points in [0, 1]
    x_grid = np.linspace(0, 1, n_points)
    
    # Start with initial density rho_0(x) = 1
    rho_vals = np.ones(n_points)

    print("Iterating the Perron-Frobenius operator...")
    for i in range(n_iter):
        # Create an interpolation function for the current density
        # Use fill_value for extrapolation, though most x_val will be small
        rho_func = interp1d(x_grid, rho_vals, kind='linear', fill_value="extrapolate")
        
        # Apply the operator to each point on the grid
        new_rho_vals = np.array([perron_frobenius_operator(rho_func, y) for y in x_grid])
        
        # Normalize the integral of the new density to 1
        integral = np.trapz(new_rho_vals, x_grid)
        if integral > 1e-9:
            rho_vals = new_rho_vals / integral
        else:
            print("Warning: Integral is close to zero. Aborting.")
            return None, None
        print(f"Iteration {i+1}/{n_iter} complete. Integral: {integral:.4f}")

    # Final normalization
    final_integral = np.trapz(rho_vals, x_grid)
    if final_integral > 1e-9:
        rho_vals /= final_integral
        
    return x_grid, rho_vals

def linear_fit(x, a, b):
    return a * x + b

if __name__ == '__main__':
    # Solve for the density
    x_grid, rho_vals = solve_invariant_density()

    if x_grid is not None:
        # Fit a line to the numerical result
        popt, _ = curve_fit(linear_fit, x_grid, rho_vals)
        a, b = popt

        # The normalized density is approximated by rho(x) = a*x + b
        print("\n---")
        print("The normalized invariant density can be approximated by a linear function.")
        print("Final Equation: rho(x) = a * x + b")
        print("With coefficients:")
        print(f"a = {a}")
        print(f"b = {b}")
        
        # Output the numbers in the final equation as requested
        print("\nNumbers in the final equation:")
        print(a)
        print(b)
        
        # Verification of the normalization for the linear fit
        integral_fit = a/2 + b
        print(f"\nThe integral of the fitted line from 0 to 1 is {integral_fit:.6f}, which should be close to 1.")
        
        # <<<>>> format for final answer: The density function itself, represented by the linear fit.
        final_answer_str = f"rho(x) = {a}*x + {b}"
        print(f"\n<<<rho(x) = {a:.4f}*x + {b:.4f}>>>")
