def find_frog_phenotype():
    """
    Identifies the correct phenotype for the Isla Colón morph of Oophaga pumilio.
    """
    # A dictionary holding the answer choices provided.
    answer_choices = {
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

    # The correct choice based on established zoological knowledge.
    correct_key = 'E'

    # Explanation based on the common name and known phenotype.
    print("The 'Isla Colón' morph of Oophaga pumilio is famously known as the 'blue jeans' poison frog.")
    print("This name directly reflects its phenotype, which is a vibrant red body with blue legs and feet.")
    print("\nMatching this description to the provided options leads to the following conclusion:")

    # Print the correct option and its full description.
    print(f"\nThe best description is Answer {correct_key}: {answer_choices[correct_key]}")

# Execute the function to find and print the answer.
find_frog_phenotype()