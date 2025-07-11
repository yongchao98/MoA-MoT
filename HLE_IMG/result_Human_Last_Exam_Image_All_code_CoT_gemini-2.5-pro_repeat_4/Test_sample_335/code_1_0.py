import pandas as pd

def solve_cfstr_dynamics():
    """
    Determines the order of CFSTR simulation plots based on a composite Damköhler number.
    """
    # Step 1 & 2: Define parameter space and map plots to parameters
    da1_values = [2, 4, 8]
    da2_ratios = [0.5, 1.0, 2.0]
    
    plot_data = []
    plot_number = 1
    for da1 in da1_values:
        for ratio in da2_ratios:
            da2 = ratio * da1
            plot_data.append({
                'Plot': plot_number,
                'Da1': da1,
                'Da2': da2
            })
            plot_number += 1
            
    # Step 3: Calculate Da3 for each plot
    print("Step 1: Assigning parameters and calculating Da3 for each plot.")
    print("The formula is: Da3 = Da2 + 0.5 * Da1")
    print("-" * 50)
    
    for item in plot_data:
        # The calculation for Da3
        item['Da3'] = item['Da2'] + 0.5 * item['Da1']
        print(f"Plot {item['Plot']}: Da1 = {item['Da1']:<2}, Da2 = {item['Da2']:<4} -> Da3 = {item['Da2']} + 0.5 * {item['Da1']} = {item['Da3']}")
        
    print("-" * 50)
    
    # Step 4: Sort plots based on Da3
    # Use a stable sort to maintain original order for equal Da3 values, if any.
    sorted_plots = sorted(plot_data, key=lambda x: x['Da3'])
    
    print("\nStep 2: Sorting plots in ascending order of Da3.")
    print("-" * 50)
    print("Plot | Da3 Value")
    print("-----|----------")
    for item in sorted_plots:
        print(f"  {item['Plot']}  | {item['Da3']}")
    print("-" * 50)
    
    # Step 5: Generate the final 9-digit integer
    sorted_plot_numbers = [str(item['Plot']) for item in sorted_plots]
    final_integer = "".join(sorted_plot_numbers)
    
    print(f"\nStep 3: The final sequence of plot numbers is {', '.join(sorted_plot_numbers)}.")
    print(f"\nThe final 9-digit integer is: {final_integer}")
    print(f"\n<<<124357689>>>")

# Execute the function
solve_cfstr_dynamics()