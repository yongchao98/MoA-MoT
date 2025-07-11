def find_frog_description():
    """
    This function identifies the correct phenotype for the 'Isla Colón' morph
    of Oophaga pumilio from a list of choices.
    """
    answer_choices = {
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

    # The key identifying features of the Isla Colón morph are a red body and blue legs.
    key_feature1 = "red body"
    key_feature2 = "blue legs"

    correct_answer_letter = None
    correct_answer_description = ""

    # Iterate through the choices to find the one that matches our key features.
    for letter, description in answer_choices.items():
        # Check if both key features are in the description text.
        if key_feature1 in description.lower() and key_feature2 in description.lower():
            correct_answer_letter = letter
            correct_answer_description = description
            break

    if correct_answer_letter:
        print(f"The best description for the 'Isla Colón' morph of Oophaga pumilio is:")
        print(f"Answer Choice {correct_answer_letter}: {correct_answer_description}")
    else:
        print("Could not find a matching description.")

find_frog_description()