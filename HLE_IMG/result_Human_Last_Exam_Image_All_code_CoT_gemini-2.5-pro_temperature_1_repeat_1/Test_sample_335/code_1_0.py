import pandas as pd

def solve_reactor_puzzle():
    """
    Determines the Da1 and Da2 values for each plot, calculates Da3,
    and sorts the plots based on Da3 to find the final 9-digit integer.
    """
    
    # Mapping rows and columns to parameter values based on analysis
    da1_map = {1: 2, 2: 4, 3: 8}  # Row -> Da1
    da2_ratio_map = {1: 0.5, 2: 1.0, 3: 2.0} # Column -> Da2/Da1 ratio

    plots_data = []
    # The plots are arranged in a 3x3 grid, numbered 1 to 9 left-to-right, top-to-bottom.
    for plot_num in range(1, 10):
        # Determine row and column from plot number
        row = (plot_num - 1) // 3 + 1
        col = (plot_num - 1) % 3 + 1
        
        # Determine Da1 and Da2
        da1 = da1_map[row]
        da2_ratio = da2_ratio_map[col]
        da2 = da2_ratio * da1
        
        # Calculate Da3 = Da2 + 1/2 * Da1
        da3 = da2 + 0.5 * da1
        
        plots_data.append({
            'plot_number': plot_num,
            'Da1': da1,
            'Da2': da2,
            'Da3': da3
        })

    # Create a DataFrame for easy viewing and sorting
    df = pd.DataFrame(plots_data)
    
    # Sort the DataFrame by the 'Da3' value in ascending order
    sorted_df = df.sort_values(by='Da3').reset_index(drop=True)
    
    print("Calculating Da3 = Da2 + 0.5 * Da1 for each plot:")
    for index, row in sorted_df.iterrows():
        # Print the calculation for each plot as requested
        print(f"Plot {int(row['plot_number'])}: Da3 = {row['Da2']} + 0.5 * {int(row['Da1'])} = {row['Da3']}")
    
    # Get the sorted plot numbers
    sorted_plot_numbers = sorted_df['plot_number'].astype(int).tolist()
    
    # Join the numbers into a single string
    final_sequence = "".join(map(str, sorted_plot_numbers))
    
    print("\nSorted plot numbers based on ascending Da3 value:")
    print(final_sequence)

solve_reactor_puzzle()