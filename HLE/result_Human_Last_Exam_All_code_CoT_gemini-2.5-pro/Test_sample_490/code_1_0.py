def find_frog_phenotype():
    """
    This function stores the descriptions of Oophaga pumilio morphs
    and identifies the correct one for the "Isla Colón" population.
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

    # The "Isla Colón" morph of Oophaga pumilio is well-documented.
    # It is famously characterized by a vibrant red body contrasted
    # by blue legs, often referred to as the "blue jeans" frog.
    correct_choice = 'E'

    # Print the correct answer and its description.
    print("The correct description for the phenotype of the 'Isla Colón' Oophaga pumilio morph is:")
    print(f"{correct_choice}: {phenotype_options[correct_choice]}")

find_frog_phenotype()