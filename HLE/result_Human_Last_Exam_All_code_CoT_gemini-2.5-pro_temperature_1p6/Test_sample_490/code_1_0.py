def find_oophaga_morph():
    """
    Analyzes descriptions of Oophaga pumilio morphs to identify the one
    matching the "Isla Colón" variant.
    """

    # The known phenotype of the "Isla Colón" morph is a bright red/orange-red body with blue legs.
    # We will score each choice based on how well it matches these key features.
    # Scoring criteria:
    # +1 point for "red" or "orange" body.
    # +1 point for "blue legs".
    # -1 point for an incorrect feature.
    print("Evaluating Oophaga pumilio morph descriptions to identify the 'Isla Colón' phenotype.")
    print("The target phenotype is a 'bright red body' with 'blue legs'.\n")

    choices = {
        'A': 'Solid black with a metallic blue sheen',
        'B': 'Orange with dark blue splotches on the legs',
        'C': 'Bright blue with black spots',
        'D': 'Green-yellow body with brown or blue legs',
        'E': 'Bright red body with blue legs and feet',
        'F': 'Solid red',
        'G': 'Purple-red body with yellow eyes',
        'H': 'Bright blue with white spots across the back',
        'I': 'Yellow with black stripes across the back'
    }

    best_choice = ''
    max_score = -999

    print("--- Scoring Analysis ---")

    for letter, description in choices.items():
        score = 0
        equation_parts = []

        # Score body color
        body_score = 0
        if 'red' in description.lower() or 'orange' in description.lower():
            body_score = 1
        elif 'blue' in description.lower() or 'yellow' in description.lower() or 'black' in description.lower() or 'green-yellow' in description.lower():
             body_score = -1 # Penalize clearly wrong body colors
        equation_parts.append(f"Body Score: {body_score}")


        # Score leg color
        leg_score = 0
        if 'blue legs' in description.lower():
            leg_score = 1
        elif 'legs' in description.lower() and 'blue' not in description.lower():
            leg_score = -1 # Penalize if legs are mentioned but are not blue
        equation_parts.append(f"Legs Score: {leg_score}")
        
        score = body_score + leg_score
        
        # Build the final equation string for printing
        equation = f"Score for Choice {letter} = {body_score} (Body) + {leg_score} (Legs)"
        print(f"{equation} = {score}. Description: '{description}'")

        if score > max_score:
            max_score = score
            best_choice = letter

    print("\n--- Conclusion ---")
    print(f"The highest score is {max_score}, which corresponds to Choice {best_choice}.")
    print(f"The description for Choice {best_choice} is: '{choices[best_choice]}'")
    print("This accurately describes the bright red body and blue legs characteristic of the Isla Colón morph.")


# Run the analysis
find_oophaga_morph()