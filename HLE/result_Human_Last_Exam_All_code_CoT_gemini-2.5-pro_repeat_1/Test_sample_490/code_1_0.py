import sys

def solve_frog_phenotype():
    """
    Identifies the correct phenotype description for the "Isla Colón"
    morph of Oophaga pumilio from a list of choices.
    """
    # The "Isla Colón" morph of Oophaga pumilio is widely known for its
    # distinct and vibrant coloration, which has earned it the nickname
    # "blue jeans" frog. This name refers to its bright red body and blue legs.
    
    known_phenotype_description = "A bright red body with blue legs."

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

    # Find the best match from the choices
    correct_key = None
    for key, description in answer_choices.items():
        if "red body" in description.lower() and "blue legs" in description.lower():
            correct_key = key
            break
            
    if correct_key:
        print(f"The 'Isla Colón' morph of Oophaga pumilio is commonly known as the 'blue jeans' frog.")
        print(f"This is because of its characteristic phenotype: {known_phenotype_description}")
        print(f"The answer choice that best describes this is:")
        print(f"{correct_key}: {answer_choices[correct_key]}")

        # Final answer output
        sys.stdout.write("\n<<<E>>>\n")
    else:
        print("Could not find a matching description.")

solve_frog_phenotype()