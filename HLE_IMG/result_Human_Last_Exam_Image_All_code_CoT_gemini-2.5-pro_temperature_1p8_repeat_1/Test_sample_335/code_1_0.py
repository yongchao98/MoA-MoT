import pandas as pd

def solve_reactor_puzzle():
    """
    Determines the order of CFSTR simulation plots based on a derived parameter Da3.
    """
    # Step 1 & 2: Define parameter space and map plots to parameters
    da1_values = [2, 4, 8]
    da2_ratios = [0.5, 1, 2]
    
    plots_data = []
    plot_number = 1
    for i in range(3):  # Corresponds to rows
        for j in range(3):  # Corresponds to columns
            da1 = da1_values[i]
            da2 = da1 * da2_ratios[j]
            
            # Step 3: Calculate Da3
            da3 = da2 + 0.5 * da1
            
            plots_data.append({
                'plot_number': plot_number,
                'Da1': da1,
                'Da2': da2,
                'Da3': da3
            })
            plot_number += 1
            
    # Step 4: Sort the plots based on Da3
    sorted_plots = sorted(plots_data, key=lambda x: x['Da3'])
    
    # Step 5: Display the sorted results and the final answer
    print("Plots arranged in ascending order of Da3 = Da2 + 1/2 * Da1:\n")
    
    final_sequence = []
    print(f"{'Plot #':<8} | {'Da3 Value':<12} | Calculation (Da2 + 0.5 * Da1)")
    print("-" * 50)
    
    for plot in sorted_plots:
        final_sequence.append(str(plot['plot_number']))
        print(f"Plot {plot['plot_number']:<5} | {plot['Da3']:<12.1f} | {plot['Da2']} + 0.5 * {plot['Da1']}")
        
    final_number = "".join(final_sequence)
    print("\nThe 9-digit integer representing the sorted plot numbers is:")
    print(final_number)

# Run the solver
solve_reactor_puzzle()