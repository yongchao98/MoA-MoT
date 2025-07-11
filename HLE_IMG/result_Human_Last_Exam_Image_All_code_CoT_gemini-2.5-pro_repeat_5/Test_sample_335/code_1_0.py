def solve_reactor_puzzle():
    """
    This script determines the order of plots based on a calculated parameter Da3.
    """
    # Step 1: Map each plot number to its (Da1, Da2) values based on visual analysis.
    # The grid is organized by increasing Da1 (top to bottom) and increasing Da2/Da1 ratio (left to right).
    plot_params = {
        1: {'da1': 2, 'da2': 1},   # Da1=2, Da2=0.5*2
        2: {'da1': 2, 'da2': 2},   # Da1=2, Da2=1.0*2
        3: {'da1': 2, 'da2': 4},   # Da1=2, Da2=2.0*2
        4: {'da1': 4, 'da2': 2},   # Da1=4, Da2=0.5*4
        5: {'da1': 4, 'da2': 4},   # Da1=4, Da2=1.0*4
        6: {'da1': 4, 'da2': 8},   # Da1=4, Da2=2.0*4
        7: {'da1': 8, 'da2': 4},   # Da1=8, Da2=0.5*8
        8: {'da1': 8, 'da2': 8},   # Da1=8, Da2=1.0*8
        9: {'da1': 8, 'da2': 16},  # Da1=8, Da2=2.0*8
    }

    # Step 2: Calculate Da3 for each plot and print the calculation.
    results = []
    print("Calculating Da3 = Da2 + 1/2 * Da1 for each plot:")
    for plot_num in sorted(plot_params.keys()):
        params = plot_params[plot_num]
        da1 = params['da1']
        da2 = params['da2']
        
        # Calculate Da3
        da3 = da2 + 0.5 * da1
        
        # Store the result for sorting
        results.append({'plot': plot_num, 'da3': da3})
        
        # Format for printing, removing .0 for whole numbers
        da3_str = f"{da3:.1f}".replace('.0', '')
        print(f"Plot {plot_num}: Da3 = {da2} + 0.5 * {da1} = {da3_str}")

    # Step 3: Sort the plots based on the calculated Da3 value.
    sorted_results = sorted(results, key=lambda x: x['da3'])

    # Step 4: Extract the sorted plot numbers.
    sorted_plot_numbers = [item['plot'] for item in sorted_results]

    # Step 5: Combine the sorted plot numbers into a single 9-digit integer.
    final_integer = "".join(map(str, sorted_plot_numbers))

    print("\nSorted plot order by Da3:")
    print(" ".join(map(str, sorted_plot_numbers)))
    
    print("\nFinal 9-digit integer:")
    print(final_integer)

solve_reactor_puzzle()