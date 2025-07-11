import textwrap

def find_frog_phenotype():
    """
    This function stores the phenotype choices for Oophaga pumilio and identifies
    the correct one for the "Isla Colón" morph.
    """

    # A dictionary mapping answer choices to their descriptions.
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

    # The "Isla Colón" morph of Oophaga pumilio is famously known for its
    # bright red or orange-red body and contrasting blue legs. This gives it
    # the common name "blue jeans" poison dart frog.
    correct_choice = 'E'

    # Print the explanation and the correct answer.
    explanation = (
        "The Oophaga pumilio morph endemic to Isla Colón in the Bocas del Toro "
        "Archipelago is widely recognized by its vibrant and distinct coloration, "
        "which has earned it the nickname 'blue jeans' frog."
    )

    print(textwrap.fill(explanation, width=80))
    print("\nBased on this, the best description is:")
    print(f"Choice {correct_choice}: {phenotypes[correct_choice]}")

find_frog_phenotype()