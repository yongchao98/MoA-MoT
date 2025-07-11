def solve_frog_phenotype():
    """
    This function identifies and prints the correct description for the
    phenotype of the 'Isla Colón' morph of the poison-dart frog Oophaga pumilio.
    """

    # A dictionary holding the answer choices provided in the problem.
    phenotype_options = {
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

    # Based on herpetological data, the Oophaga pumilio from Isla Colón
    # is famous for its "blue jeans" appearance.
    # This corresponds to a vividly colored body and blue legs.
    correct_choice_letter = 'E'

    print("Analyzing the phenotypes of Oophaga pumilio morphs...")
    print("The 'Isla Colón' morph is well-documented.")
    print("Its distinctive coloration gives it the common name 'blue jeans' frog.")
    print("\nThe correct description is:")
    
    # Printing the 'equation' of the answer by showing the letter and its value
    print(f"Choice {correct_choice_letter}: {phenotype_options[correct_choice_letter]}")

solve_frog_phenotype()