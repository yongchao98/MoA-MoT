import pandas as pd

def solve_reactor_puzzle():
    """
    Solves the CFSTR plot identification puzzle by:
    1. Defining the mapping from plot number to Damköhler numbers (Da1, Da2) based on visual analysis.
    2. Calculating the sorting metric Da3 = Da2 + 0.5 * Da1 for each plot.
    3. Sorting the plots in ascending order of their Da3 values.
    4. Printing the final 9-digit integer result.
    """

    # Step 1: Map plots to (Da1, Da2) pairs based on visual analysis of reactor dynamics.
    # The keys are plot numbers (1-9), and values are tuples of (Da1, Da2).
    # This mapping is derived by matching system behavior (stable vs. oscillatory,
    # relative reactant consumption) to the corresponding parameter values.
    plot_to_da_map = {
        1: (2, 2),    # Small oscillations, Y ≈ Z -> Da1=Da2, lowest oscillatory reactivity
        2: (4, 8),    # Oscillatory, Z < Y -> Da2 > Da1
        3: (8, 16),   # Most violent oscillations, Z < Y -> Da2 > Da1, highest reactivity
        4: (4, 4),    # Oscillatory, Y ≈ Z -> Da1=Da2
        5: (8, 8),    # Oscillatory, Y ≈ Z -> Da1=Da2, more reactive than plot 4
        6: (2, 4),    # Stable, Z < Y -> Da2 > Da1
        7: (8, 4),    # Oscillatory, Y < Z -> Da1 > Da2
        8: (4, 2),    # Stable, Y < Z -> Da1 > Da2
        9: (2, 1),    # Stable, least reactive -> lowest Da values
    }

    plot_data = []
    # Step 2: Calculate Da3 for each plot
    for plot_num, (da1, da2) in plot_to_da_map.items():
        da3 = da2 + 0.5 * da1
        plot_data.append({'plot': plot_num, 'Da1': da1, 'Da2': da2, 'Da3': da3})

    # Create a DataFrame for easy sorting and display
    df = pd.DataFrame(plot_data)

    # Step 3: Sort the DataFrame by the 'Da3' column
    df_sorted = df.sort_values(by='Da3').reset_index(drop=True)

    print("Analysis Results: Plot assignments and calculated Da3 values")
    print(df_sorted.to_string())
    print("\n" + "="*50 + "\n")

    # Step 4: Generate the final 9-digit number
    # The final equation is the sequence of plot numbers arranged in ascending order of Da3.
    # Here we output each number in that final sequence.
    final_sequence = df_sorted['plot'].tolist()
    final_number_str = "".join(map(str, final_sequence))
    
    print("The final sequence of plot numbers, sorted by Da3, is composed of the following digits:")
    for number in final_sequence:
        print(f"Digit: {number}")
        
    print("\nFinal 9-digit integer:")
    print(final_number_str)


solve_reactor_puzzle()