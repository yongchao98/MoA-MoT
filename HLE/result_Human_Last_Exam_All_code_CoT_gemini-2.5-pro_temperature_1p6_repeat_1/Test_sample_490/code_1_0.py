def find_frog_phenotype():
    """
    Identifies and prints the correct description for the Oophaga pumilio 'Isla Colón' morph.
    """
    # A dictionary mapping answer choices to their descriptions.
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

    # The Oophaga pumilio morph from Isla Colón is famously known for its
    # bright red or orange-red body and contrasting blue legs, earning it
    # the nickname "blue jeans" poison frog.
    correct_choice = 'E'

    # Print the final answer components as requested.
    print("The best description for the phenotype of the Oophaga pumilio 'Isla Colón' morph is:")
    print(f"Answer Choice: {correct_choice}")
    print(f"Description: {phenotype_choices[correct_choice]}")

find_frog_phenotype()