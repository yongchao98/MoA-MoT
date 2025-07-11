def solve_tapestry_mystery():
    """
    This function simulates the microbe's effect on the tapestry
    to determine the mutated enzyme and original patch color for an orange result.
    """
    
    # Let's lay out all possible scenarios
    original_colors = ["yellow", "blue"]
    mutated_enzymes = ["A", "B", "C", "D"]
    
    solution = None
    
    # We will loop through every possibility
    for color in original_colors:
        for enzyme in mutated_enzymes:
            
            # Start with the original pigment
            pigments = {color}
            
            # Define which enzymes are functional
            has_A = enzyme != "A"
            has_B = enzyme != "B"
            has_C = enzyme != "C"
            has_D = enzyme != "D"
            
            # --- Simulate the metabolic process ---
            
            # Pathway 1: Blue degradation
            if "blue" in pigments and has_D:
                pigments.remove("blue")
                pigments.add("yellow")

            # Pathway 2: Yellow degradation
            final_color = ""
            if "yellow" in pigments and has_A:
                # If Enzyme B is missing, the Red intermediate accumulates.
                # The problem states the result is orange, meaning both red and yellow are present.
                # This implies the Yellow -> Red conversion is incomplete when B is missing.
                if not has_B:
                    # This is our key condition for orange
                    final_color = "orange"
                else:
                    # Otherwise, the process continues
                    pigments.remove("yellow")
                    pigments.add("red_intermediate")
            
            if "red_intermediate" in pigments and has_B:
                pigments.remove("red_intermediate")
                pigments.add("blue_intermediate")
                
            if "blue_intermediate" in pigments and has_C:
                pigments.remove("blue_intermediate")
                # Final product is colorless
                pigments = {"colorless"}
            
            # Determine final color if not already orange
            if not final_color:
                final_color = list(pigments)[0] if pigments else "unknown"

            # Check if this scenario results in the observed orange color
            if final_color == "orange":
                # We assume the simplest path is the intended solution
                # The path from yellow to orange is more direct than blue -> yellow -> orange
                if color == "yellow":
                    solution = f"{enzyme}-{color}"
    
    print(f"The mutated enzyme and original color must be: {solution}")

solve_tapestry_mystery()
<<<B-yellow>>>