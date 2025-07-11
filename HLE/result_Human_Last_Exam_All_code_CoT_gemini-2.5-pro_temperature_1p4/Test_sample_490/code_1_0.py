def find_frog_phenotype():
    """
    Identifies and explains the correct phenotype for the Oophaga pumilio 'Isla Colón' morph.
    """
    
    # A dictionary mapping the answer choice letters to their descriptions.
    phenotype_options = {
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

    # The correct answer based on established herpetological knowledge.
    correct_answer_key = 'E'

    # Retrieve the text of the correct answer.
    correct_description = phenotype_options[correct_answer_key]

    print("Analysis of the Oophaga pumilio 'Isla Colón' Morph:")
    print("=" * 50)
    print("The frog Oophaga pumilio, commonly known as the strawberry poison-dart frog, exhibits a wide range of color morphs across its habitat.")
    print("The specific morph endemic to Isla Colón in the Bocas del Toro Archipelago of Panamá is one of the most famous.")
    print("\nThis morph is characterized by its striking and distinct coloration:")
    print(" - Body: A vibrant, bright red.")
    print(" - Legs and Feet: A contrasting deep blue or purplish-blue.")
    print("\nThis pattern is often nicknamed the 'blue jeans' morph.")
    print("=" * 50)
    
    # Final conclusion based on the analysis.
    print(f"Comparing this information with the given choices, the best description is:")
    print(f"\nFinal Answer choice ({correct_answer_key}): {correct_description}")

# Execute the function to provide the answer.
find_frog_phenotype()