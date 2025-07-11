def find_frog_phenotype():
    """
    This function identifies the correct phenotype for the
    Oophaga pumilio "Isla Colón" morph from a list of choices.
    """
    # Description of choices provided to the user.
    choices = {
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

    # The Oophaga pumilio morph from Isla Colón in the Bocas del Toro Archipelago
    # is famously known as the "blue jeans" frog.
    # Its phenotype is a bright red or orange-red body with dark blue or purplish-blue legs.
    correct_choice_key = 'E'
    
    print("Analyzing the provided phenotypes for the Oophaga pumilio 'Isla Colón' morph...")
    print(f"The correct description is: '{choices[correct_choice_key]}'")
    print(f"This corresponds to answer choice: {correct_choice_key}")

find_frog_phenotype()