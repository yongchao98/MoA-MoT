def solve_tapestry_puzzle():
    """
    This function systematically evaluates the scenarios to find the cause of the orange patch.
    """
    possible_original_colors = ["yellow", "blue"]
    possible_missing_enzymes = ["A", "B", "C", "D"]
    solution_found = False

    for original_color in possible_original_colors:
        for missing_enzyme in possible_missing_enzymes:
            final_color = ""

            # Case 1: Original patch was yellow
            # Pathway: yellow --(A)--> red --(B)--> blue_intermediate --(C)--> colorless
            if original_color == "yellow":
                if missing_enzyme == 'A':
                    final_color = "yellow"  # First step blocked
                elif missing_enzyme == 'B':
                    # Red intermediate accumulates. Mix with original yellow makes orange.
                    final_color = "orange"
                elif missing_enzyme == 'C':
                    final_color = "blue"    # Blue intermediate accumulates
                elif missing_enzyme == 'D':
                    final_color = "yellow"  # Enzyme D is not in this pathway

            # Case 2: Original patch was blue
            # Pathway: blue --(D)--> yellow
            elif original_color == "blue":
                if missing_enzyme == 'D':
                    final_color = "blue"   # Pathway blocked
                else: # Missing A, B, or C
                    final_color = "yellow" # Enzyme D is present and converts blue to yellow

            # Check if this scenario creates an orange patch
            if final_color == "orange":
                # Print the components of the answer as requested
                print(f"The mutated enzyme is: {missing_enzyme}")
                print(f"The original colour of the patch was: {original_color}")

                # Print the final answer in the required format
                final_answer = f"{missing_enzyme}-{original_color}"
                print(f"Final formatted answer: {final_answer}")
                solution_found = True
                break
        if solution_found:
            break

solve_tapestry_puzzle()