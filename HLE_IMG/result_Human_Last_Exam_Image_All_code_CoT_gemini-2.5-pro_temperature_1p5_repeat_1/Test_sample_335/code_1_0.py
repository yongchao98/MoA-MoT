import collections

def solve_reactor_puzzle():
    """
    Determines the order of reactor simulation plots based on a calculated parameter Da3.
    """
    
    # Step 1: Define the parameter values.
    # Da1 values are associated with rows (top to bottom).
    # Da2/Da1 ratios are associated with columns (left to right).
    da1_values = [2, 4, 8]
    da2_ratios = [0.5, 1.0, 2.0]

    # Step 2: Map each plot number to its corresponding (Da1, Da2) pair
    # and calculate Da3 for each.
    # The grid is 3x3, numbered 1-3 (top), 4-6 (middle), 7-9 (bottom).
    plots_data = []
    plot_number = 1
    for da1 in da1_values:
        for ratio in da2_ratios:
            da2 = ratio * da1
            
            # The formula for the sorting key.
            # Da3 = Da2 + 1/2 * Da1
            da3 = da2 + 0.5 * da1
            
            plots_data.append({
                'plot_num': plot_number,
                'da1': da1,
                'da2': da2,
                'da3': da3
            })
            plot_number += 1

    # Print the parameters and the calculation for each plot as requested.
    print("Calculating Da3 for each plot based on its position in the grid:")
    print("-" * 55)
    print(f"{'Plot':<6} | {'Da1':<5} | {'Da2':<6} | {'Calculation of Da3':<25}")
    print("-" * 55)
    for plot in plots_data:
        # Here we output each number in the final equation.
        equation_str = f"{plot['da2']:.1f} + 0.5 * {plot['da1']} = {plot['da3']:.1f}"
        print(f"{plot['plot_num']:<6} | {plot['da1']:<5} | {plot['da2']:<6.1f} | {equation_str:<25}")
    print("-" * 55)
    
    # Step 3: Sort the plots based on their Da3 value in ascending order.
    sorted_plots = sorted(plots_data, key=lambda p: p['da3'])
    
    # Step 4: Extract the plot numbers from the sorted list and create the final answer.
    sorted_plot_numbers = [str(p['plot_num']) for p in sorted_plots]
    final_integer_string = "".join(sorted_plot_numbers)
    
    print("\nPlots sorted by ascending Da3 value:")
    sorted_sequence = " -> ".join(sorted_plot_numbers)
    print(sorted_sequence)
    
    print("\nFinal 9-digit integer:")
    print(final_integer_string)

# Execute the function to solve the puzzle.
solve_reactor_puzzle()