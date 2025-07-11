def solve_frog_phenotype_question():
    """
    This function analyzes the provided options and identifies the correct
    phenotype for the 'Isla Colón' morph of the Oophaga pumilio frog.
    """

    question = 'Which of the following best describes the phenotype of the "Isla Colón" morph of Oophaga pumilio?'

    # A dictionary to hold the answer choices
    options = {
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

    # The "Isla Colón" morph is famously characterized by its bright red body
    # and blue or dark-blue legs, often called the "blue jeans" frog.
    correct_choice_key = 'E'

    print("Question:")
    print(question)
    print("\nEvaluating Options:")

    for key, description in options.items():
        print(f"  {key}: {description}")

    print("\n------------------------------------------------\n")
    print("Conclusion:")
    print("The 'Isla Colón' morph of Oophaga pumilio is one of the most recognized variants.")
    print("It is known for its distinct and vibrant coloration which serves as a warning to predators.")
    print("\nThe correct description is:")
    print(f"Answer: [{correct_choice_key}] {options[correct_choice_key]}")


solve_frog_phenotype_question()