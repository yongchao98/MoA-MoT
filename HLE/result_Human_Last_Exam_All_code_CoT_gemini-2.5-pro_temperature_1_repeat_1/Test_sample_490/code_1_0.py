def find_frog_phenotype():
    """
    This function identifies and prints the correct phenotype for the
    'Isla Colón' morph of Oophaga pumilio from a list of choices.
    """
    # A dictionary mapping the answer choices to their descriptions.
    phenotype_choices = {
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

    # The "Isla Colón" morph of Oophaga pumilio is known for its
    # iconic, solid bright red coloration.
    correct_answer_key = 'F'

    # Retrieve the description of the correct answer.
    correct_description = phenotype_choices[correct_answer_key]

    # Print the final answer.
    print(f"The correct phenotype description is:")
    print(f"Choice {correct_answer_key}: {correct_description}")

find_frog_phenotype()