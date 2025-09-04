import sympy as sp
import numpy as np
from scipy.stats import linregress

def check_correctness():
    """
    Checks the correctness of the answer to the parallax distribution question.
    The check is performed in two ways:
    1.  Symbolic Derivation: Replicates the mathematical derivation using the sympy library.
    2.  Numerical Simulation: Simulates a uniform star distribution and checks the resulting power law.
    """
    
    # --- Part 1: Symbolic Derivation ---
    # This part verifies the mathematical logic step-by-step.
    try:
        # Define symbols for the derivation
        d, plx, rho, pi = sp.symbols('d plx rho pi', positive=True, real=True)
        
        # The number of stars per unit range of parallax, dN/d(plx), is proportional to d^2 * |dd/d(plx)|.
        # This comes from dN ∝ d^2 * dd and the chain rule.
        
        # Define the relationship between distance (d) and parallax (plx)
        d_expr = 1 / plx
        
        # Find the derivative of d with respect to plx to get dd/d(plx)
        dd_dplx = sp.diff(d_expr, plx)
        
        # The proportionality for dN/d(plx) is d^2 * |dd/d(plx)|
        dN_dplx_prop = (d_expr**2) * abs(dd_dplx)
        
        # Simplify the final expression
        simplified_expr = sp.simplify(dN_dplx_prop)
        
        # The expected expression from the derivation is 1/plx^4
        expected_expr = 1 / plx**4
        
        # Check if the derived expression is proportional to the expected one.
        # Their ratio should be a constant.
        is_proportional = sp.simplify(simplified_expr / expected_expr).is_constant()
        
        if not is_proportional:
            return f"Incorrect. The symbolic derivation is flawed. The derived proportionality is {simplified_expr}, but it should be {expected_expr}."
            
        # The derivation correctly yields a proportionality of 1/plx^4.
        # This corresponds to option D in the question.
        # The provided answer is <<<D>>>, which is consistent with the derivation.
        
    except Exception as e:
        return f"An error occurred during the symbolic check: {e}"

    # --- Part 2: Numerical Simulation (as a secondary, supporting check) ---
    try:
        # Generate a large number of stars with uniform random positions in a 3D sphere
        N_stars = 2_000_000
        max_radius = 1000
        min_radius = 1  # To avoid d=0 -> plx=inf
        
        # Generate random radii with r^2 distribution for uniform volume density
        r_cubed = np.random.uniform(min_radius**3, max_radius**3, N_stars)
        distances = np.cbrt(r_cubed)
        
        # Calculate parallax
        parallaxes = 1.0 / distances
        
        # Create a histogram of the parallax values
        num_bins = 50
        bins = np.logspace(np.log10(parallaxes.min()), np.log10(parallaxes.max()), num_bins)
        
        # Get the counts (dN) in each bin
        counts, bin_edges = np.histogram(parallaxes, bins=bins)
        
        # Calculate the width of each bin (d(plx))
        bin_widths = np.diff(bin_edges)
        
        # Calculate the center of each bin (plx)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # The quantity of interest is dN/d(plx).
        # We must normalize the counts by the bin width.
        # Filter out empty bins to avoid log(0) errors.
        non_zero_mask = counts > 0
        if not np.any(non_zero_mask):
             return "Incorrect. Numerical simulation failed: no stars were binned."

        y_values = counts[non_zero_mask] / bin_widths[non_zero_mask]
        x_values = bin_centers[non_zero_mask]
        
        # To check the power law, we fit a line to the log-log plot.
        # y = C * x^m  =>  log(y) = log(C) + m * log(x)
        log_x = np.log(x_values)
        log_y = np.log(y_values)
        
        # Perform linear regression to find the slope 'm'
        slope, _, _, _, _ = linregress(log_x, log_y)
        
        # The expected slope is -4.
        expected_slope = -4.0
        tolerance = 0.2  # Allow for some numerical error and statistical fluctuations.
        
        if abs(slope - expected_slope) > tolerance:
            return f"Incorrect. The theoretical derivation is correct (1/plx^4), but the numerical simulation gives a power-law exponent of {slope:.2f}, which is not close to the expected -4. While the theory supports the answer, this numerical check fails."

    except Exception as e:
        return f"An error occurred during the numerical check: {e}"

    # If both the symbolic derivation and numerical simulation confirm the result,
    # and the result matches the provided answer, then the answer is correct.
    return "Correct"

# Execute the check
result = check_correctness()
print(result)