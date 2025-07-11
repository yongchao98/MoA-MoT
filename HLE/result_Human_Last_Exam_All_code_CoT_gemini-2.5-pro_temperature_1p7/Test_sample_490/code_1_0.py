import sys

def find_frog_phenotype():
    """
    This function identifies the correct phenotype for the Oophaga pumilio "Isla Colón" morph.
    """
    question = "Which of the following best describes the phenotype of the 'Isla Colón' morph of Oophaga pumilio?"

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

    # The "Isla Colón" morph of Oophaga pumilio is widely documented.
    # Its distinct and famous coloration is a vibrant red body with blue legs.
    # This combination matches choice E.
    correct_answer_key = 'E'

    print(f"Question: {question}")
    print("\nAnalysis:")
    print("The Oophaga pumilio morph endemic to Isla Colón in the Bocas del Toro Archipelago is known for its aposematic (warning) coloration.")
    print("This specific population famously exhibits a bright red body contrasted with blue legs and feet.")
    print("\nConclusion:")
    print(f"The correct choice is '{correct_answer_key}', which is: '{answer_choices[correct_answer_key]}'.")


find_frog_phenotype()