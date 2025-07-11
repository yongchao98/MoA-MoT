def solve_tapestry_puzzle():
    """
    Solves the microbe tapestry puzzle by simulating each possible scenario.
    The goal is to find which scenario results in an 'orange' color.
    """

    print("Analyzing the metabolic pathways for each mutation scenario...\n")

    target_color = 'orange'
    solution_found = False
    final_answer_enzyme = ""
    final_answer_color = ""

    original_colors = ['yellow', 'blue']
    mutations = ['A', 'B', 'C', 'D']

    # Define color mixing rules. Orange is a mix of yellow and red.
    color_map = {
        frozenset(['yellow']): 'yellow',
        frozenset(['red']): 'red',
        frozenset(['blue_pigment']): 'blue',
        frozenset(['blue_intermediate']): 'blue',
        frozenset(['colourless']): 'colourless',
        frozenset(['yellow', 'red']): 'orange',
    }

    # Iterate through all possibilities
    for start_color in original_colors:
        for missing_enzyme in mutations:
            
            # This set will store the pigments that accumulate at the end
            final_pigments = set()

            if start_color == 'yellow':
                # Pathway: yellow ->(A) red ->(B) blue_intermediate ->(C) colourless
                if missing_enzyme == 'A':
                    # Process is blocked at the start, yellow pigment remains.
                    final_pigments.add('yellow')
                elif missing_enzyme == 'B':
                    # Yellow is converted to red, which then accumulates.
                    final_pigments.add('red')
                elif missing_enzyme == 'C':
                    # The pathway proceeds until the blue intermediate, which accumulates.
                    final_pigments.add('blue_intermediate')
                elif missing_enzyme == 'D':
                    # Enzyme D is not in this pathway. The process completes to the end.
                    final_pigments.add('colourless')
            
            elif start_color == 'blue':
                # Combined pathway: blue ->(D) yellow ->(A) red ->(B) blue_intermediate ->(C) colourless
                if missing_enzyme == 'D':
                    # Process is blocked at the start, blue pigment remains.
                    final_pigments.add('blue_pigment')
                elif missing_enzyme == 'A':
                    # Blue is converted to yellow, which then accumulates as it cannot proceed.
                    final_pigments.add('yellow')
                elif missing_enzyme == 'B':
                    # Blue is converted to yellow. Yellow is converted to red. Red accumulates.
                    # An orange color implies a mix of yellow and red are present.
                    # This happens when yellow is being produced while red is accumulating.
                    final_pigments.add('yellow')
                    final_pigments.add('red')
                elif missing_enzyme == 'C':
                    # The entire pathway runs until the final step. Blue intermediate accumulates.
                    final_pigments.add('blue_intermediate')

            # Determine final color from the set of pigments
            final_color_name = color_map.get(frozenset(final_pigments), 'unknown')

            # Print the step-by-step thinking for this scenario
            print(f"Scenario:")
            print(f"  - Original Patch Color: {start_color.capitalize()}")
            print(f"  - Lacking Enzyme: {missing_enzyme}")
            pigment_str = ", ".join(sorted(list(final_pigments))) if final_pigments else "none"
            print(f"  - Resulting Pigments: {{{pigment_str}}}")
            print(f"  - Final Color: {final_color_name.capitalize()}")
            
            if final_color_name == target_color:
                solution_found = True
                final_answer_enzyme = missing_enzyme
                final_answer_color = start_color
                print("  - This matches the target color (Orange)!")

            print("-" * 25)

    if solution_found:
        final_answer_str = f"{final_answer_enzyme}-{final_answer_color}"
        print(f"\nConclusion:")
        print(f"The only scenario resulting in an orange color is an originally '{final_answer_color}' patch with a microbe lacking enzyme '{final_answer_enzyme}'.")
        print(f"\nThe equation is: Missing Enzyme = {final_answer_enzyme}, Original Color = {final_answer_color}")
        print(f"<<<{final_answer_str}>>>")
    else:
        print("\nCould not find a scenario that results in an orange patch based on the rules.")

solve_tapestry_puzzle()