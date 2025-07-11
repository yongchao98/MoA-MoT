import collections

def solve_roller_puzzle():
    """
    This function solves the roller drive puzzle by matching configurations to displacement plots.
    The logic is based on analyzing the number of lobes and the amplitude of variation.
    """

    # Step 1: Characterize each roller configuration (1-8).
    # 'lobes': Number of lobes on the driving (green) roller, determining oscillation frequency.
    # 'variation': A qualitative score (1=low, 2=medium, 3=high) for the expected amplitude of speed variation
    #              based on the roller shapes.
    configs = {
        1: {'lobes': 6, 'variation': 3},  # High variation: Deep driver lobes, pointy driven lobes.
        2: {'lobes': 6, 'variation': 1},  # Low variation: Shallow driver lobes dampen the effect of the pointy driven roller.
        3: {'lobes': 2, 'variation': 3},  # High variation: Pronounced lobes and concave/convex shapes.
        4: {'lobes': 4, 'variation': 3},  # High variation: Very pointy and sharp shapes on both rollers.
        5: {'lobes': 6, 'variation': 2},  # Medium variation: Deep driver lobes but a more rounded driven roller than config 1.
        6: {'lobes': 4, 'variation': 1},  # Low variation: Both rollers are smooth with shallow lobes/shape.
        7: {'lobes': 5, 'variation': 2},  # Medium variation: Moderately defined shapes.
        8: {'lobes': 3, 'variation': 2},  # Medium variation: Moderately defined shapes.
    }

    # Step 2: Characterize each displacement plot (A-H).
    # 'oscillations': Number of "wiggles" observed in the plot.
    # 'amplitude': A qualitative score for the observed variation in the curve's slope.
    plots = {
        'A': {'oscillations': 2, 'amplitude': 3},
        'B': {'oscillations': 6, 'amplitude': 3},
        'C': {'oscillations': 5, 'amplitude': 2},
        'D': {'oscillations': 4, 'amplitude': 3},
        'E': {'oscillations': 6, 'amplitude': 2},
        'F': {'oscillations': 6, 'amplitude': 1},
        'G': {'oscillations': 4, 'amplitude': 1},
        'H': {'oscillations': 3, 'amplitude': 2},
    }

    # Step 3: Match plots to configurations based on the defined characteristics.
    pairings = {}
    available_configs = list(configs.keys())

    # Iterate through plots in alphabetical order to ensure a consistent result.
    for plot_letter in sorted(plots.keys()):
        plot_data = plots[plot_letter]
        
        # Find a configuration that has the same number of lobes and variation score.
        for config_num in available_configs:
            config_data = configs[config_num]
            
            if (config_data['lobes'] == plot_data['oscillations'] and
                config_data['variation'] == plot_data['amplitude']):
                
                pairings[plot_letter] = config_num
                available_configs.remove(config_num)  # Ensure each config is used only once.
                break

    # Step 4: Print the results as requested.
    print("The determined pairings (the 'equation') are:")
    final_sequence_list = []
    for plot_letter in sorted(pairings.keys()):
        config_number = pairings[plot_letter]
        # This line outputs each number in the final "equation" as requested.
        print(f"Plot {plot_letter} = Configuration {config_number}")
        final_sequence_list.append(str(config_number))

    final_sequence = "".join(final_sequence_list)
    print("\nThe final sequence of configuration numbers for plots A-H is:")
    print(final_sequence)

solve_roller_puzzle()