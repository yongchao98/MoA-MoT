import pandas as pd

def solve_reactor_puzzle():
    """
    Solves the CFSTR plot puzzle by mapping parameters, calculating Da3,
    sorting the plots, and printing the results.
    """
    # Step 1: Establish the mapping between plots and (Da1, Da2) parameters.
    # Based on the analysis of the plots' dynamic behavior.
    plot_params = {
        1: {'Da1': 2, 'Da2': 1},
        2: {'Da1': 2, 'Da2': 2},
        3: {'Da1': 2, 'Da2': 4},
        4: {'Da1': 4, 'Da2': 2},
        5: {'Da1': 4, 'Da2': 4},
        6: {'Da1': 4, 'Da2': 8},
        7: {'Da1': 8, 'Da2': 4},
        8: {'Da1': 8, 'Da2': 8},
        9: {'Da1': 8, 'Da2': 16}
    }

    # Step 2: Calculate Da3 for each plot and store the results.
    results = []
    for plot_num, params in plot_params.items():
        da1 = params['Da1']
        da2 = params['Da2']
        da3 = da2 + 0.5 * da1
        results.append({
            'plot': plot_num,
            'da1': da1,
            'da2': da2,
            'da3': da3
        })

    # Step 3: Sort the plots based on the calculated Da3 value.
    sorted_results = sorted(results, key=lambda x: x['da3'])

    # Step 4: Print the ordered calculations.
    print("The plots are sorted by the value of Da3 = Da2 + 0.5 * Da1.")
    print("Here are the calculations in ascending order of Da3:")
    for res in sorted_results:
        print(f"Plot {res['plot']}: Da3 = {res['da2']} + 0.5 * {res['da1']} = {res['da3']}")

    # Step 5: Construct and print the final 9-digit integer.
    final_integer = "".join([str(res['plot']) for res in sorted_results])
    print(f"\nThe resulting 9-digit integer is:")
    print(final_integer)

solve_reactor_puzzle()