def solve_tapestry_puzzle():
    """
    Simulates the effect of each possible single-enzyme mutation
    on patches of yellow or blue to find the cause of the orange patch.
    """
    
    # Possible original colours and missing enzymes
    original_colours = ["yellow", "blue"]
    missing_enzymes = ["A", "B", "C", "D"]
    
    print("Analyzing all possible scenarios...\n")
    
    # This dictionary will hold the final colour for each scenario
    results = {}
    
    for colour in original_colours:
        for enzyme in missing_enzymes:
            scenario = f"Original: {colour}, Missing: {enzyme}"
            
            # Check for presence of each enzyme
            has_A = (enzyme != "A")
            has_B = (enzyme != "B")
            has_C = (enzyme != "C")
            has_D = (enzyme != "D")
            
            final_colour = ""
            
            if colour == "yellow":
                # Pathway: Yellow -> Red -> Blue_Intermediate -> Colourless
                if not has_A:
                    final_colour = "yellow" # Blocked at the start
                elif not has_B:
                    final_colour = "red" # Y -> R, then blocked
                elif not has_C:
                    final_colour = "blue_intermediate" # Y -> R -> Bi, then blocked
                else: # has A, B, C. Missing D has no effect.
                    final_colour = "colourless" # Pathway runs to completion
            
            elif colour == "blue":
                # Pathway: Blue -> Yellow -> Red -> Blue_Intermediate -> Colourless
                if not has_D:
                    final_colour = "blue" # Blocked at the start
                else:
                    # Blue is converted to Yellow by D. Now analyze fate of Yellow.
                    if not has_A:
                        final_colour = "yellow" # B -> Y, then blocked
                    elif not has_B:
                        # B -> Y (via D), Y -> R (via A)
                        # This is the special case. The production of yellow from the
                        # blue source and the conversion of yellow to red happen
                        # simultaneously, resulting in a mix of yellow and red.
                        final_colour = "orange" 
                    elif not has_C:
                        final_colour = "blue_intermediate" # B -> Y -> R -> Bi, then blocked
                    else: # has A, B, C, D
                        final_colour = "colourless" # Entire pathway runs
            
            results[scenario] = final_colour
            print(f"- {scenario:<25} -> Final Colour: {final_colour}")
            
    print("\nFinding the scenario that results in 'orange'...")
    
    # Find the correct answer
    final_enzyme = None
    final_original_colour = None
    for scenario, result_colour in results.items():
        if result_colour == "orange":
            # Parsing the key string like "Original: blue, Missing: B"
            parts = scenario.split(',')
            final_original_colour = parts[0].split(': ')[1]
            final_enzyme = parts[1].split(': ')[1]
            break
            
    if final_enzyme and final_original_colour:
        answer = f"{final_enzyme}-{final_original_colour}"
        print(f"\nConclusion: The patch is orange when the missing enzyme is {final_enzyme} and the original colour was {final_original_colour}.")
        print("Final Answer:")
        print(answer)
        # The final answer in the required format
        print(f"\n<<<{answer}>>>")
    else:
        print("Could not find a scenario that results in an orange patch.")

solve_tapestry_puzzle()