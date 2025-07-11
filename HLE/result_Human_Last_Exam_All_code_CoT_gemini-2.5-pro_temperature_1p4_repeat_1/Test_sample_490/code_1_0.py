def find_frog_phenotype():
    """
    This function identifies the correct description for the Oophaga pumilio
    morph from Isla Colón, Panamá.
    """
    # A dictionary mapping the answer choices to their descriptions.
    phenotype_options = {
        'A': "Solid black with a metallic blue sheen",
        'B': "Orange with dark blue splotches on the legs",
        'C': "Bright blue with black spots",
        'D': "Green-yellow body with brown or blue legs",
        'E': "Bright red body with blue legs and feet",
        'F': "Solid red",
        'G': "Purple-red body with yellow eyes",
        'H': "Bright blue with white spots across the back",
        'I': "Yellow with black stripes across the back"
    }

    # Based on herpetological knowledge, the Isla Colón morph of Oophaga pumilio
    # is famous for its bright red body and blue limbs.
    correct_answer_key = 'E'

    print("Analyzing the phenotypes for Oophaga pumilio, Isla Colón morph...")
    print("-" * 30)
    print(f"The correct description is found in option {correct_answer_key}.")
    print(f"Description: {phenotype_options[correct_answer_key]}")
    print("-" * 30)
    print("This morph is commonly known for its 'strawberry' red body and 'blue jeans' legs, making this the most accurate choice.")

# Execute the function to print the answer.
find_frog_phenotype()