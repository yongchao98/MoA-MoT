def solve_color_puzzle():
    """
    This script determines the color of a concentrated Pseudomonas aeruginosa sample
    by explaining the pigments it produces.
    """

    # 1. Define the pigments and their colors
    pigment_colors = {
        "Pyocyanin": "Blue",
        "Pyoverdine": "Yellow-Green"
    }

    # 2. The combination of these pigments gives the final color
    final_color = "Blue-Green"

    # 3. Print the explanation
    print("Step 1: Pseudomonas aeruginosa produces key pigments.")
    print(f" - The pigment Pyocyanin is {pigment_colors['Pyocyanin']}.")
    print(f" - The pigment Pyoverdine is {pigment_colors['Pyoverdine']}.")
    print("\nStep 2: The color of the bacterial sample is due to the combination of these pigments.")

    # 4. Print the final "equation" as requested, showing each component
    print("\nFinal Color Equation:")
    component1 = pigment_colors['Pyocyanin']
    component2 = pigment_colors['Pyoverdine']
    result = final_color
    print(f"'{component1}' + '{component2}' => '{result}'")

    print(f"\nConclusion: A dense, concentrated sample of Pseudomonas aeruginosa is characteristically {final_color}.")

solve_color_puzzle()