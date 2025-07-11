def solve_reactor_puzzle():
    """
    This script determines the Da1 and Da2 values for each of the 9 plots,
    calculates a composite value Da3, and sorts the plots based on Da3
    to produce a final 9-digit integer.
    """
    # Step 1: Define the mapping of plot numbers to (Da1, Da2) pairs.
    # This is based on the analysis that Da1 increases across columns (2, 4, 8)
    # and the ratio Da2/Da1 increases down rows (0.5, 1, 2).
    plot_params = {
        # Row 1: Da2 = 0.5 * Da1
        1: (2, 1),
        2: (4, 2),
        3: (8, 4),
        # Row 2: Da2 = 1.0 * Da1
        4: (2, 2),
        5: (4, 4),
        6: (8, 8),
        # Row 3: Da2 = 2.0 * Da1
        7: (2, 4),
        8: (4, 8),
        9: (8, 16),
    }

    # Step 2: Calculate Da3 = Da2 + 0.5 * Da1 for each plot.
    # The results will be stored in a list of tuples: (Da3_value, plot_number)
    print("Calculating Da3 = Da2 + 0.5 * Da1 for each plot:")
    plot_da3_values = []
    for plot_num in sorted(plot_params.keys()):
        da1, da2 = plot_params[plot_num]
        da3 = da2 + 0.5 * da1
        plot_da3_values.append((da3, plot_num))
        # Output the equation for each plot as requested
        print(f"Plot {plot_num}: Da3 = {da2} + 0.5 * {da1} = {da3}")

    # Step 3: Sort the plots in ascending order based on the Da3 value.
    plot_da3_values.sort()

    # Step 4: Extract the sorted plot numbers and create the final 9-digit integer.
    print("\nPlots sorted by Da3 in ascending order:")
    print([f"Plot {p[1]} (Da3={p[0]})" for p in plot_da3_values])
    
    sorted_plot_numbers = [item[1] for item in plot_da3_values]
    final_integer_str = "".join(map(str, sorted_plot_numbers))

    print(f"\nThe final sequence of plot numbers is: {final_integer_str}")

solve_reactor_puzzle()