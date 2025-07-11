def find_morph_phenotype():
    """
    This function identifies the correct phenotype for the 'Isla Colón'
    morph of Oophaga pumilio from a list of choices.
    """
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

    # Based on biological data, the 'Isla Colón' morph is known for a red body and blue legs.
    # We will select the option that matches this description.
    correct_answer_key = 'E'

    correct_description = phenotype_options[correct_answer_key]

    print("The question asks for the phenotype of the 'Isla Colón' morph of Oophaga pumilio.")
    print("This morph is famously known for its 'blue jeans' appearance.")
    print("\nSearching the options for the best match...")
    print(f"\nFinal Answer Found:")
    print(f"Choice: {correct_answer_key}")
    print(f"Description: {correct_description}")

find_morph_phenotype()