def find_frog_phenotype():
    """
    This function stores the possible phenotypes for the Oophaga pumilio morph
    and identifies the correct one for the "Isla Colón" morph.
    """
    phenotypes = {
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

    # The "Isla Colón" morph of Oophaga pumilio is well-documented
    # as having a bright red body with blue legs.
    correct_answer_key = 'E'

    print("The following options were considered:")
    for key, value in phenotypes.items():
        print(f"{key}: {value}")

    print("\nBased on biological data, the correct description is:")
    print(f"Answer {correct_answer_key}: {phenotypes[correct_answer_key]}")

find_frog_phenotype()