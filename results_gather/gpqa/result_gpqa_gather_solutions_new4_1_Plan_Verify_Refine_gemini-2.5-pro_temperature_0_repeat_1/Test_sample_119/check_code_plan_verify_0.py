import numpy as np
from scipy.stats import linregress

def check_answer():
    """
    This function checks the correctness of the answer by running a numerical simulation.
    
    The question asks how the number of stars per unit parallax range (dN/d(plx)) changes with parallax (plx),
    assuming a uniform 3D distribution of stars.

    The theoretical derivation, as correctly shown in the provided answer's analysis, leads to:
    dN/d(plx) ∝ 1/plx^4

    This means a log-log plot of dN/d(plx) vs. plx should be a straight line with a slope of -4.
    This simulation will:
    1. Generate a large number of stars uniformly in a 3D spherical shell.
    2. Calculate the parallax for each star (plx = 1/distance).
    3. Create a histogram of the parallaxes to approximate dN/d(plx).
    4. Perform a linear regression on the log-log plot of the histogram data.
    5. Check if the resulting slope is close to -4, which would validate the answer A (∝ 1/plx^4).
    """
    try:
        # 1. Define simulation parameters
        num_stars = 2_000_000
        # Define a spherical shell to avoid parallaxes that are too large (d->0) or too small (d->inf)
        min_distance = 1.0
        max_distance = 1000.0

        # 2. Generate stars with a uniform 3D distribution within the shell
        # To get a uniform distribution in a sphere, the probability density of finding a star
        # at radius r is proportional to r^2. We use inverse transform sampling on the volume.
        # The cumulative distribution function for r is proportional to r^3.
        # r = (R1^3 + u * (R2^3 - R1^3))^(1/3), where u is uniform in [0,1].
        u = np.random.rand(num_stars)
        distances = np.cbrt(min_distance**3 + u * (max_distance**3 - min_distance**3))

        # 3. Calculate parallaxes (plx = 1/d)
        parallaxes = 1.0 / distances

        # 4. Create a histogram to find dN/d(plx)
        # The range of parallaxes is from 1/max_distance to 1/min_distance.
        num_bins = 100
        counts, bin_edges = np.histogram(parallaxes, bins=num_bins)
        
        # Calculate bin centers and the (constant) bin width
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_width = bin_edges[1] - bin_edges[0]

        # Calculate dN/d(plx) = counts / bin_width
        # Filter out bins with zero counts to avoid issues with log(0)
        non_zero_mask = counts > 0
        if np.sum(non_zero_mask) < 10: # Need enough points for a reliable fit
             return "Simulation failed: Not enough data points for a reliable fit."

        counts_nz = counts[non_zero_mask]
        bin_centers_nz = bin_centers[non_zero_mask]
        
        n_plx = counts_nz / bin_width

        # 5. Perform a log-log fit to check the power law
        # Expected relationship: n(plx) = C * plx^-4
        # In log space: log(n(plx)) = log(C) - 4 * log(plx)
        # This is a linear relationship y = c + m*x, where y=log(n(plx)), x=log(plx), and m=-4.
        log_x = np.log(bin_centers_nz)
        log_y = np.log(n_plx)

        # Perform linear regression
        slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)

        # 6. Check the result against the provided answer
        # The provided answer is A, which corresponds to 1/plx^4, so the expected slope is -4.
        expected_slope = -4.0
        tolerance = 0.15 # Allow for some numerical error in the simulation

        # Check if the slope is close to the expected value and the fit is good (high R-squared)
        if abs(slope - expected_slope) < tolerance and r_value**2 > 0.99:
            # The simulation confirms the theoretical result 1/plx^4, which is option A.
            # The provided answer is also A. Therefore, the answer is correct.
            return "Correct"
        else:
            # The simulation result does not match the expected value for answer A.
            reason = (f"Incorrect. The provided answer is A, which corresponds to the relationship dN/d(plx) ∝ 1/plx^4. "
                      f"This implies that a log-log plot of the number of stars per unit parallax vs. parallax should have a slope of -4. "
                      f"However, a numerical simulation of the scenario yields a slope of {slope:.4f} with an R-squared value of {r_value**2:.4f}. "
                      f"This numerical result does not support the provided answer.")
            return reason

    except Exception as e:
        return f"An error occurred during the checking code: {e}"

# To run the check, you would execute the following line:
# print(check_answer())