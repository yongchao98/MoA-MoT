import math

def solve_reactor_puzzle():
    """
    Determines the Da1 and Da2 values for each plot, calculates Da3,
    and sorts the plots based on Da3 to generate a 9-digit integer.
    """
    # Parameter space
    da1_values = [2, 4, 8]  # Corresponds to rows 1, 2, 3
    da2_ratios = [0.5, 1.0, 2.0]  # Corresponds to columns 1, 2, 3

    plot_data = []

    print("Calculating Da3 = Da2 + 1/2 * Da1 for each plot:")
    print("-" * 50)

    # Iterate through plots 1 to 9 to calculate Da3 for each
    for plot_number in range(1, 10):
        # Determine row and column from plot number (0-indexed)
        row_index = (plot_number - 1) // 3
        col_index = (plot_number - 1) % 3

        # Assign Da1 and Da2 based on the grid layout
        da1 = da1_values[row_index]
        da2 = da1 * da2_ratios[col_index]

        # Calculate Da3
        da3 = da2 + 0.5 * da1

        # Store the results
        plot_data.append({'plot': plot_number, 'da1': da1, 'da2': da2, 'da3': da3})

        # Print the equation for each plot
        print(f"Plot {plot_number}: Da3 = {da2:.1f} + 0.5 * {da1:.1f} = {da3:.1f}")

    # Sort the list of dictionaries by the 'da3' value
    sorted_plots = sorted(plot_data, key=lambda x: x['da3'])

    # Extract the plot numbers from the sorted list
    sorted_plot_numbers = [item['plot'] for item in sorted_plots]

    # Combine the sorted plot numbers into a single string
    final_integer_string = "".join(map(str, sorted_plot_numbers))

    print("-" * 50)
    print("Plots sorted by ascending Da3 value:", sorted_plot_numbers)
    print("\nFinal 9-digit integer:")
    print(final_integer_string)
    
    # The final answer format required by the system
    print(f"\n<<<{final_integer_string}>>>")

solve_reactor_puzzle()