def solve_reactor_puzzle():
    """
    This function solves the CFSTR plot puzzle by:
    1. Mapping plots 1-9 to their corresponding (Da1, Da2) parameters based on a grid hypothesis.
    2. Calculating Da3 for each plot using the formula Da3 = Da2 + 0.5 * Da1.
    3. Sorting the plots based on the calculated Da3 value.
    4. Printing the final 9-digit integer sequence.
    """
    # Define the parameter space
    da1_values = [2, 4, 8]  # Corresponds to rows 1, 2, 3
    da2_ratios = [0.5, 1, 2] # Corresponds to columns 1, 2, 3

    plot_data = []

    print("Calculating Da3 for each plot based on the hypothesized parameter mapping:")
    # Iterate through the 3x3 grid to assign parameters to each plot
    for i in range(3):  # Row index
        for j in range(3):  # Column index
            plot_number = i * 3 + j + 1
            da1 = float(da1_values[i])
            da2 = float(da2_ratios[j] * da1)

            # Calculate Da3
            da3 = da2 + 0.5 * da1

            plot_data.append({
                'plot_number': plot_number,
                'Da1': da1,
                'Da2': da2,
                'Da3': da3
            })
            # Output the equation for each plot
            print(f"Plot {plot_number} (Da1={int(da1)}, Da2={int(da2)}): Da3 = {da2} + 1/2 * {da1} = {da3}")

    # Sort the plots by the Da3 value in ascending order
    sorted_plots = sorted(plot_data, key=lambda x: x['Da3'])

    # Extract the plot numbers from the sorted list
    sorted_plot_numbers = [str(p['plot_number']) for p in sorted_plots]

    # Join the numbers to form the final 9-digit integer string
    final_sequence = "".join(sorted_plot_numbers)

    print("\nThe sequence of plot numbers arranged in ascending order of Da3 is:")
    print(final_sequence)

solve_reactor_puzzle()