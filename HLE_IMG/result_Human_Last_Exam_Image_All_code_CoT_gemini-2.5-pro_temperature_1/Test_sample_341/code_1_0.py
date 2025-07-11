import math

def solve_statistical_mechanics_plots():
    """
    This function solves the problem by identifying the plots and calculating R_max.
    """

    # Step 1 & 2: Identify g(r) plots for each system in the specified order {SS, SR, R, HS, TW}
    # SS (Plot 1): Repulsive shoulder creates a dip between r=1 and r=1.5.
    # SR (Plot 3): Stickiness creates a delta-function spike at r=1.
    # R (Plot 5): Repulsive ramp creates an increasing slope from r=1 to r=1.5.
    # HS (Plot 7): Hard core creates a simple jump at r=1.
    # TW (Plot 9): Attractive well creates a high peak at r=1, followed by a decay.
    g_r_indices = [1, 3, 5, 7, 9]

    # Step 3 & 4: Identify S(k) plots and the unique system.
    # S(0) analysis: Attraction increases S(0), repulsion decreases it vs HS.
    # S_SR (Plot 6): Highest S(0) due to strong stickiness.
    # S_TW (Plot 4): High S(0) due to attractive well.
    # S_HS (Plot 8): Reference S(0) value.
    # S_R (Plot 2): Low S(0) due to repulsion. The g(r) for R (Plot 5) is smoother than for SS (Plot 1),
    # which corresponds to faster decaying S(k) oscillations as seen in Plot 2.
    # This leaves SS as the unique system with no S(k) plot.
    s_k_indices = [0, 6, 2, 8, 4] # Order: {SS, SR, R, HS, TW}

    # The final 10-element list of plot indices.
    plot_indices = g_r_indices + s_k_indices

    # Step 5: Calculate R_max for the unique system (SS, plot 1).
    # R_g(r) = g(r+1)/g(r) for r = {1.5, 2.5, 3.5, ...} (in units of sigma).
    # r=1.5: R_g(1.5) = g(2.5)/g(1.5). From plot 1, g(1.5) is high (~2.0), g(2.5) is the first trough (~0.9). Ratio is < 1.
    # r=2.5: R_g(2.5) = g(3.5)/g(2.5). g(3.5) is the second peak (~1.05), g(2.5) is the first trough (~0.9).
    # The ratio is ~1.05/0.9 = 1.1666..., which is 7/6. This is the maximum value.
    # r=3.5: R_g(3.5) = g(4.5)/g(3.5). This ratio is again < 1.
    R_max = 7.0 / 6.0

    # Step 6: Assemble and print the final result.
    final_values = plot_indices + [R_max]

    # Print the result in the format {v1, v2, ..., v11}
    print("{", end="")
    for i in range(len(final_values) - 1):
        value = final_values[i]
        if isinstance(value, int):
            print(f"{value}, ", end="")
        else:
            print(f"{value}, ", end="")
    print(f"{final_values[-1]}}}", end="")

solve_statistical_mechanics_plots()