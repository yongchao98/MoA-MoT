def find_pumilio_morph_phenotype():
    """
    Analyzes choices to identify the phenotype of the Oophaga pumilio 'Isla Colón' morph.
    """
    # Step 1: Define the provided answer choices.
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

    # Step 2: Define the key identifying features of the Isla Colón morph.
    # This morph is famously known for its "blue jeans" appearance.
    key_feature_body = "red body"
    key_feature_legs = "blue legs"

    print("Analyzing phenotypes to identify the 'Isla Colón' morph...")
    print("=========================================================")
    print(f"The target phenotype has two key features: a '{key_feature_body}' and '{key_feature_legs}'.")

    # Step 3 & 4: Evaluate each choice and formulate the equation.
    best_choice = None
    max_score = 0

    for letter, description in choices.items():
        current_score = 0
        if key_feature_body in description.lower():
            current_score += 1
        if key_feature_legs in description.lower():
            current_score += 1

        if current_score > max_score:
            max_score = current_score
            best_choice = letter

    # Formulate and print the scoring equation for the best match.
    score_for_red_body = 1
    score_for_blue_legs = 1
    total_score = score_for_red_body + score_for_blue_legs

    print("\nScoring Equation for a perfect match:")
    print(f"Match for '{key_feature_body}' + Match for '{key_feature_legs}' = Total Score")
    print(f"       {score_for_red_body}          +          {score_for_blue_legs}          =    {total_score}")
    print("=========================================================")

    # Step 5: Output the result.
    print(f"\nThe choice that perfectly matches this scoring equation is Choice {best_choice}.")
    print(f"Final Answer: The description for the Isla Colón morph is '{choices[best_choice]}'.")


find_pumilio_morph_phenotype()