import sys

def solve_tapestry_puzzle():
    """
    Simulates the microbe metabolism to find the cause of the orange patch.
    """
    original_colours = ["yellow", "blue"]
    missing_enzymes = ["A", "B", "C", "D"]
    solution_found = False

    print("Analyzing the metabolic pathways to find the cause of the orange patch...\n")

    # Iterate through each possible original colour
    for colour in original_colours:
        # Iterate through each possible single enzyme mutation
        for enzyme in missing_enzymes:
            if solution_found:
                break

            print(f"--- Analyzing Case: Original Colour = {colour.capitalize()}, Missing Enzyme = {enzyme} ---")
            
            final_pigments = set()
            explanation = ""

            if colour == "yellow":
                # Pathway: yellow -> red -> blue_intermediate -> colourless
                if enzyme == "A":
                    # Enzyme A is missing, so yellow pigment is not degraded.
                    final_pigments.add("yellow")
                    explanation = "Enzyme A is missing. The yellow pigment cannot be converted to the red intermediate. The patch remains yellow."
                else:
                    # Enzyme A is present, yellow is converted to red.
                    if enzyme == "B":
                        # Enzyme B is missing, red intermediate accumulates.
                        # The reaction yellow -> red is not instantaneous, resulting in a mix.
                        final_pigments.add("yellow")
                        final_pigments.add("red")
                        explanation = "Enzyme B is missing. Yellow pigment is converted to the red intermediate, but the pathway stops there. The mix of remaining yellow pigment and red intermediate appears orange."
                    else:
                        # Enzyme B is present, red is converted to blue_intermediate.
                        if enzyme == "C":
                            # Enzyme C is missing, blue_intermediate accumulates.
                            final_pigments.add("blue_intermediate")
                            explanation = "Enzyme C is missing. The pathway proceeds to the blue intermediate and stops. The patch turns blue."
                        else:
                            # Enzymes A, B, C are present. Pathway completes.
                            final_pigments.add("colourless")
                            explanation = "The pathway is fully functional (Enzyme D mutation is irrelevant). The yellow pigment is degraded to a colourless product."
            
            elif colour == "blue":
                # Pathway: blue -> yellow
                if enzyme == "D":
                    # Enzyme D is missing, blue pigment is not degraded.
                    final_pigments.add("blue")
                    explanation = "Enzyme D is missing. The blue pigment cannot be degraded. The patch remains blue."
                else:
                    # Enzyme D is present, blue is converted to yellow.
                    # Missing enzymes A, B, or C do not affect this pathway.
                    final_pigments.add("yellow")
                    explanation = f"Enzyme D is functional. The blue pigment is converted to a yellow product. The patch turns yellow."

            # Determine the final colour from the set of pigments
            final_colour = ""
            if "red" in final_pigments and "yellow" in final_pigments:
                final_colour = "orange"
            elif "red" in final_pigments:
                final_colour = "red"
            elif "yellow" in final_pigments:
                final_colour = "yellow"
            elif "blue_intermediate" in final_pigments or "blue" in final_pigments:
                final_colour = "blue"
            elif "colourless" in final_pigments:
                final_colour = "colourless"

            print(f"Explanation: {explanation}")
            print(f"Resulting Pigments: {final_pigments}")
            print(f"Final Colour: {final_colour}\n")

            if final_colour == "orange":
                print("="*40)
                print("SOLUTION FOUND!")
                print(f"The patch is orange because the microbe is missing Enzyme '{enzyme}'.")
                print(f"This caused the original '{colour}' pigment to be partially converted to a red intermediate,")
                print("resulting in a mix of red and yellow pigments.")
                print("="*40)
                print(f"\nFinal Answer: {enzyme}-{colour}")
                solution_found = True
                # Use sys.exit() to stop the script after finding the solution
                # In a real script, you might return the value instead.
                # For this interactive environment, we'll let the loops finish.
                # sys.exit() 

if __name__ == "__main__":
    solve_tapestry_puzzle()