def solve_tapestry_mystery():
    """
    This function logically deduces the mutated enzyme and original tapestry color.
    It simulates each possible scenario and identifies the one resulting in an orange patch.
    """
    
    mutations = ['A', 'B', 'C', 'D']
    original_colors = ['yellow', 'blue']
    solution = None

    print("Analyzing the metabolic pathways to find the cause of the orange patch...")
    print("-" * 60)

    for color in original_colors:
        for enzyme in mutations:
            print(f"Scenario: Original color is {color.upper()}, mutated enzyme is {enzyme.upper()}")
            
            final_color = ""
            explanation = ""

            if color == 'yellow':
                # Pathway: yellow ->(A)-> red ->(B)-> blue_intermediate ->(C)-> colorless
                if enzyme == 'A':
                    final_color = 'yellow'
                    explanation = "Enzyme A is required for the first step. With a mutation in A, the yellow pigment is not metabolized. Result: The patch remains YELLOW."
                elif enzyme == 'B':
                    final_color = 'orange'
                    explanation = "Enzyme A converts yellow pigment to a red intermediate. Enzyme B is required for the next step. With a mutation in B, the red intermediate accumulates. The mix of remaining original yellow pigment and the new red intermediate makes the patch ORANGE."
                elif enzyme == 'C':
                    final_color = 'blue_intermediate'
                    explanation = "The pathway proceeds to the blue intermediate (yellow -> red -> blue_intermediate), but stops there due to the mutation in C. Result: The patch turns BLUE (the color of the intermediate)."
                elif enzyme == 'D':
                    final_color = 'colorless'
                    explanation = "Enzyme D is not in this pathway. The pathway proceeds to completion (yellow -> ... -> colorless). Result: The patch becomes COLORLESS."
            
            elif color == 'blue':
                # Pathway: blue ->(D)-> yellow
                if enzyme in ['A', 'B', 'C']:
                    final_color = 'yellow'
                    explanation = f"Enzyme {enzyme} is not in the blue pigment pathway. The pathway (blue -> yellow) is unaffected and proceeds to completion. Result: The patch turns YELLOW."
                elif enzyme == 'D':
                    final_color = 'blue'
                    explanation = "Enzyme D is required to metabolize the blue pigment. With a mutation in D, the blue pigment is unaffected. Result: The patch remains BLUE."
            
            print(f"   - Outcome: {explanation}")
            if final_color == 'orange':
                solution = f"{enzyme}-{color}"
            print("-" * 60)

    if solution:
        print("\nConclusion:")
        print("The only scenario that results in an orange patch is a mutation in Enzyme B on an originally yellow patch.")
        print("The final answer is:")
        print(f"<<<{solution}>>>")
    else:
        print("\nConclusion: Could not find a scenario that results in an orange patch.")

solve_tapestry_mystery()