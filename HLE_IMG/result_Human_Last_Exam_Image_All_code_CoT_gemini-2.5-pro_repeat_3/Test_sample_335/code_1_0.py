import pandas as pd

def solve_reactor_dynamics():
    """
    Analyzes reactor simulation plots to determine parameters and sort by a derived value.
    """
    # Step 1 & 2: Based on qualitative analysis of the plots, we map each plot
    # number to its corresponding (Da1, Da2) parameter pair.
    # This mapping is the result of analyzing stability, oscillation intensity,
    # and relative reactant consumption in each plot.
    # - Plot 1: Mildest oscillations -> Lowest Da values (2,1).
    # - Plot 9: Most stable, highest conversion -> Highest Da values (8,16).
    # - Plot 8: Stable, high Da1, lower Da2 -> (8,4) -> Y (blue) consumed more.
    # - Plot 6: Stable, high Da2, lower Da1 -> (4,8) -> Z (green) consumed more.
    # - Plots with Y~Z (blue~green lines): (2,2) -> Plot 2; (4,4) -> Plot 4; (8,8) -> Plot 3.
    # - Plot with Y<Z (blue<green): Da1>Da2 -> (4,2) -> Plot 5.
    # - Plot with Z<Y (green<blue): Da2>Da1 -> (2,4) -> Plot 7.
    
    plot_assignments = {
        1: {'Da1': 2, 'Da2': 1},
        2: {'Da1': 2, 'Da2': 2},
        3: {'Da1': 8, 'Da2': 8},
        4: {'Da1': 4, 'Da2': 4},
        5: {'Da1': 4, 'Da2': 2},
        6: {'Da1': 4, 'Da2': 8},
        7: {'Da1': 2, 'Da2': 4},
        8: {'Da1': 8, 'Da2': 4},
        9: {'Da1': 8, 'Da2': 16}
    }

    # Step 3: Calculate Da3 for each plot and store the results.
    plot_data = []
    for plot_num, params in plot_assignments.items():
        da1 = params['Da1']
        da2 = params['Da2']
        # The equation to calculate Da3
        da3 = da2 + 0.5 * da1
        plot_data.append({'plot': plot_num, 'Da1': da1, 'Da2': da2, 'Da3': da3})

    # Use pandas for easy sorting and display.
    df = pd.DataFrame(plot_data)

    # Step 4: Sort the DataFrame by the 'Da3' column in ascending order.
    df_sorted = df.sort_values(by='Da3').reset_index(drop=True)

    print("Derivation of the sorted order based on Da3 = Da2 + 1/2 * Da1:\n")
    
    # Print the calculation for each plot in the sorted order
    for index, row in df_sorted.iterrows():
        print(f"Plot {int(row['plot'])}: Da3 = {row['Da2']:>2.0f} + 0.5 * {row['Da1']:>2.0f} = {row['Da3']:.1f}")

    # Step 5: Generate the final 9-digit integer.
    sorted_plots_list = df_sorted['plot'].astype(int).tolist()
    final_answer = "".join(map(str, sorted_plots_list))
    
    print("\nPlot numbers arranged in ascending order of Da3:")
    print(final_answer)

solve_reactor_dynamics()