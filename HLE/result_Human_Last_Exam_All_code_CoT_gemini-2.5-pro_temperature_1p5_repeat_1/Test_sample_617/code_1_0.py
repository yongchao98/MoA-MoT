def get_pseudomonas_color():
    """
    Determines the color of a dense Pseudomonas aeruginosa sample based on its pigments.
    """
    # Pseudomonas aeruginosa produces two primary pigments.
    pigment1_name = "Pyocyanin"
    pigment1_color = "Blue"

    pigment2_name = "Pyoverdine"
    pigment2_color = "Yellow-Green"

    # The combination of these pigments gives the bacteria their characteristic color.
    # The prompt requires an 'equation', so we'll represent the color mixing as one.
    final_color = "Blue-Green"

    print(f"Bacterium: Pseudomonas aeruginosa")
    print(f"Pigment 1: {pigment1_name} ({pigment1_color})")
    print(f"Pigment 2: {pigment2_name} ({pigment2_color})")
    print("-" * 20)
    print("Color 'Equation':")
    # Fulfilling the requirement to output each part of the 'equation'.
    print(f"The color '{pigment1_color}' from Pyocyanin + the color '{pigment2_color.split('-')[1]}' from Pyoverdine => '{final_color}'")
    print("-" * 20)
    print(f"A dense, washed, and concentrated sample of Pseudomonas aeruginosa is colored by its pigments, resulting in a distinct '{final_color}' appearance.")
    print("Therefore, the correct answer is C.")

get_pseudomonas_color()