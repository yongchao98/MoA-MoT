import pandas as pd

def solve_reactor_dynamics():
    """
    This function determines the Da1 and Da2 values for each plot, calculates Da3,
    and orders the plots based on ascending Da3 values.
    """
    # Step 1: Define the parameter space for Da1 and the ratios for Da2.
    da1_values = [2, 4, 8]
    da2_ratios = [0.5, 1, 2]

    # Step 2: Establish the mapping from parameter pairs to plot numbers.
    # This mapping is derived from a qualitative analysis of the plots' dynamic
    # behavior (stable vs. oscillatory, amplitude, frequency) and matching it
    # to the quantitative reactivity measure Da3.
    # Format: (Da1, Da2): plot_number
    parameter_to_plot_mapping = {
        # Stable states, ordered by increasing conversion (reactivity)
        (2, 1): 8,   # Da3=2.0: Lowest reactivity, highest reactant concentration.
        (2, 2): 6,   # Da3=3.0: Medium reactivity among stable states.
        (4, 2): 9,   # Da3=4.0: Highest reactivity among stable states.
        
        # Oscillatory states, ordered by increasing intensity
        (2, 4): 1,   # Da3=5.0: Onset of instability, small amplitude oscillations.
        (4, 4): 4,   # Da3=6.0: Medium amplitude oscillations.
        (8, 4): 2,   # Da3=8.0: More intense medium amplitude oscillations.
        (4, 8): 3,   # Da3=10.0: Large amplitude, low-frequency oscillations.
        (8, 8): 7,   # Da3=12.0: More intense large amplitude oscillations.
        (8, 16): 5,  # Da3=20.0: Most intense, lowest frequency oscillations.
    }

    # Step 3: Calculate Da3 for each parameter set and associate with its plot number.
    plot_data = []
    for da1 in da1_values:
        for ratio in da2_ratios:
            da2 = da1 * ratio
            
            # Use integer keys for the dictionary if Da2 is an integer
            da2_key = int(da2) if da2.is_integer() else da2
            
            plot_number = parameter_to_plot_mapping.get((da1, da2_key))
            
            if plot_number is not None:
                # The problem defines Da3 = Da2 + 1/2*Da1
                da3 = da2 + 0.5 * da1
                plot_data.append({
                    'Plot Number': plot_number,
                    'Da1': da1,
                    'Da2': da2_key,
                    'Da3': da3
                })

    # Step 4: Sort the plots based on the calculated Da3 value in ascending order.
    sorted_plots = sorted(plot_data, key=lambda x: x['Da3'])

    # Step 5: Display the results and formulate the final answer.
    print("The plots are ordered by the ascending value of Da3 = Da2 + 0.5 * Da1.")
    print("-" * 65)
    print(f"{'Plot Number':<12} | {'Da1':<5} | {'Da2':<5} | {'Calculation':<25} | {'Da3':<5}")
    print("-" * 65)

    final_sequence = []
    for p in sorted_plots:
        final_sequence.append(str(p['Plot Number']))
        calculation_str = f"{p['Da2']} + 0.5 * {p['Da1']}"
        print(f"{p['Plot Number']:<12} | {p['Da1']:<5} | {p['Da2']:<5} | {calculation_str:<25} | {p['Da3']:<5.1f}")
    
    final_integer = "".join(final_sequence)
    
    print("-" * 65)
    print(f"\nThe sequence of plot numbers arranged by ascending Da3 is: {', '.join(final_sequence)}")
    print(f"The final 9-digit integer is: {final_integer}")

solve_reactor_dynamics()