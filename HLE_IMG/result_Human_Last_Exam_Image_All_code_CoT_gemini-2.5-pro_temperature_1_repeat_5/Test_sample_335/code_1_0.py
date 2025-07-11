import pandas as pd

def solve_reactor_puzzle():
    """
    Solves the CFSTR puzzle by mapping parameters to plots, calculating Da3,
    and sorting the plots accordingly.
    """
    # Step 1 & 2: Define parameters and map them to plots
    da1_values = [2, 4, 8]
    da2_ratios = [0.5, 1, 2]
    
    plot_data = []
    plot_number = 1
    for i, da1 in enumerate(da1_values):
        for j, ratio in enumerate(da2_ratios):
            da2 = ratio * da1
            plot_data.append({'plot': plot_number, 'Da1': da1, 'Da2': da2})
            plot_number += 1
            
    # Step 3: Calculate Da3 for each plot
    print("Calculating Da3 = Da2 + 0.5 * Da1 for each plot:")
    for data in plot_data:
        # The prompt requires showing each number in the equation
        da1 = data['Da1']
        da2 = data['Da2']
        half_da1 = 0.5 * da1
        da3 = da2 + half_da1
        data['Da3'] = da3
        print(f"Plot {data['plot']}: Da3 = {da2:.1f} + 0.5 * {da1:.1f} = {da2:.1f} + {half_da1:.1f} = {da3:.1f}")

    # Step 4: Sort the plots based on the Da3 value
    sorted_plots = sorted(plot_data, key=lambda x: x['Da3'])
    
    print("\nPlots sorted by ascending Da3 value:")
    df = pd.DataFrame(sorted_plots)
    print(df.to_string(index=False))

    # Step 5: Generate the final 9-digit integer
    final_sequence = "".join([str(p['plot']) for p in sorted_plots])
    
    print("\nFinal sequence of plot numbers arranged by ascending Da3:")
    print(final_sequence)

solve_reactor_puzzle()