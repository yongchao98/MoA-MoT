def solve_tapestry_mystery():
    """
    This function systematically analyzes the possible scenarios to determine
    the missing enzyme and original tapestry colour.
    """

    # Define the color mixing rules
    # Orange is a mix of Red and Yellow
    target_colors = {"red", "yellow"}

    # Define the metabolic pathways
    # Pathway 1: yellow -> red -> blue_intermediate -> colourless (Enzymes A, B, C)
    # Pathway 2: blue -> yellow (Enzyme D)

    possible_mutations = ['A', 'B', 'C', 'D']
    original_colours = ['yellow', 'blue']
    solution_found = False

    print("Analyzing scenarios to find the cause of the orange patch...")
    print("An orange patch means both 'red' and 'yellow' pigments must be present.\n")

    for start_colour in original_colours:
        if solution_found:
            break
        print(f"--- Testing starting colour: {start_colour.capitalize()} ---")
        for missing_enzyme in possible_mutations:
            
            # Simulate the reactions
            present_pigments = set()

            if start_colour == 'yellow':
                if missing_enzyme == 'A':
                    present_pigments = {'yellow'}
                    explanation = "Pathway stops. Patch remains yellow."
                elif missing_enzyme == 'B':
                    present_pigments = {'red'}
                    explanation = "Yellow is converted to Red, pathway stops. Patch turns red."
                elif missing_enzyme == 'C':
                    # Assuming full conversion at each step
                    present_pigments = {'blue_intermediate'}
                    explanation = "Yellow is converted to Blue Intermediate, pathway stops. Patch turns blue."
                else: # Missing D
                    present_pigments = {'colourless'}
                    explanation = "Enzyme D is not in this pathway. Yellow is fully degraded to a colourless product."
            
            elif start_colour == 'blue':
                if missing_enzyme == 'D':
                    present_pigments = {'blue'}
                    explanation = "Pathway to degrade Blue pigment is blocked. Patch remains blue."
                else:
                    # Blue is converted to Yellow by enzyme D
                    # Now we see what happens to this new yellow pigment
                    if missing_enzyme == 'A':
                        present_pigments = {'yellow'}
                        explanation = "Blue is converted to Yellow, but Pathway 1 is blocked. Patch turns yellow."
                    elif missing_enzyme == 'B':
                        # This is the key scenario
                        present_pigments = {'yellow', 'red'}
                        explanation = "Blue is converted to Yellow (by D), which is then converted to Red (by A). The pathway stops. The mix of Yellow and Red creates an orange color."
                    elif missing_enzyme == 'C':
                        present_pigments = {'blue_intermediate'}
                        explanation = "Blue is converted to Yellow, which is then fully converted to Blue Intermediate. Patch turns blue."
            
            # Check for a solution
            if present_pigments == target_colors:
                print(f"[*] Mutation: Missing Enzyme {missing_enzyme}")
                print(f"    - Result: {explanation}")
                print(f"\nConclusion: The microbe is missing enzyme {missing_enzyme} and the original patch was {start_colour}.")
                
                final_answer = f"{missing_enzyme}-{start_colour}"
                print("\nFinal Answer Equation:")
                print("Original Patch = Blue")
                print("Enzyme D converts Blue -> Yellow")
                print("Enzyme A converts Yellow -> Red")
                print("Enzyme B is MISSING, so Red cannot be converted further.")
                print("Final state = Yellow + Red = Orange")
                print(f"Therefore, the answer is: {final_answer}")
                print(f"\n<<<B-blue>>>")

                solution_found = True
                break
            else:
                 print(f"[-] Mutation: Missing Enzyme {missing_enzyme}")
                 print(f"    - Result: {explanation} Final color(s): {', '.join(list(present_pigments)) or 'none'}. No match.")
        print("")


solve_tapestry_mystery()