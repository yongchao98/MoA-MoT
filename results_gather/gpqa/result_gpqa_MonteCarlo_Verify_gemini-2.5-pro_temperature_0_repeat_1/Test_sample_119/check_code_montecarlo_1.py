import numpy as np
from scipy.stats import linregress

def check_star_parallax_distribution():
    """
    Checks the correctness of the answer D) ~ 1/plx^4 by simulating a
    uniform star distribution and fitting the resulting parallax data.
    """
    # The answer 'D' implies a power-law relationship with an exponent of -4.
    expected_exponent = -4.0

    try:
        # --- Monte Carlo Simulation ---
        # 1. Setup simulation parameters
        num_stars = 4_000_000
        max_distance = 1000.0  # Arbitrary distance units

        # 2. Generate star distances for a uniform 3D volume density.
        # To achieve this, r^3 must be sampled from a uniform distribution,
        # which means r is the cube root of a uniform random variable.
        u = np.random.rand(num_stars)
        distances = max_distance * np.cbrt(u)

        # 3. Calculate parallaxes (plx = 1/d)
        parallaxes = 1.0 / distances

        # 4. Create a log-binned histogram to find the number density n(plx).
        # Logarithmic bins are ideal for analyzing power-law distributions.
        num_bins = 50
        # Define a robust histogram range to avoid bins with poor statistics
        plx_min = 1.0 / max_distance
        plx_max = np.percentile(parallaxes, 99.5)

        bins = np.logspace(np.log10(plx_min), np.log10(plx_max), num_bins + 1)
        counts, bin_edges = np.histogram(parallaxes, bins=bins)

        # 5. Calculate the density n(plx) = dN/d(plx) ≈ ΔN/Δplx
        bin_widths = np.diff(bin_edges)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

        # Filter out bins with low counts for a more stable fit
        stable_mask = counts > 10
        if np.sum(stable_mask) < 5:
            return "Simulation failed: Not enough data points for a reliable fit."

        density = counts[stable_mask] / bin_widths[stable_mask]
        plx_values = bin_centers[stable_mask]

        # 6. Fit a power law by performing linear regression on the log-log data.
        # If density = A * plx^B, then log(density) = log(A) + B * log(plx).
        # The slope of this line is the desired exponent B.
        log_plx = np.log10(plx_values)
        log_density = np.log10(density)
        
        slope, _, r_value, _, _ = linregress(log_plx, log_density)
        simulated_exponent = slope

        # --- Verification ---
        # Check if the simulated exponent is close to the theoretical exponent (-4)
        # and if the fit quality (R-squared) is high.
        is_good_fit = (r_value**2) > 0.998
        tolerance = 0.1

        if abs(simulated_exponent - expected_exponent) < tolerance and is_good_fit:
            return "Correct"
        else:
            return (
                f"The answer 'D' (~ 1/plx^4) is incorrect based on this check.\n"
                f"Reasoning:\n"
                f"1. The answer implies a power-law exponent of {expected_exponent}.\n"
                f"2. The simulation yielded an exponent of {simulated_exponent:.4f}.\n"
                f"3. The difference of {abs(simulated_exponent - expected_exponent):.4f} is outside the tolerance of {tolerance}.\n"
                f"4. The theoretical derivation strongly supports an exponent of -4, so this discrepancy likely points to simulation noise, but the check fails based on the set criteria."
            )

    except Exception as e:
        return f"An error occurred during the checking process: {e}"

# Run the check
result = check_star_parallax_distribution()
print(result)