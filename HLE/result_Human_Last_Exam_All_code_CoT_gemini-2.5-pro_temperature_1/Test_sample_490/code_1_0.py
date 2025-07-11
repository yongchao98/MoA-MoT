import textwrap

def find_phenotype_description():
    """
    This function identifies the correct description for the Oophaga pumilio 'Isla Colón' morph.
    """
    answer_choices = {
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

    # Based on biological and herpetological sources, the "Isla Colón" morph of
    # Oophaga pumilio is famously known for its bright red to reddish-orange body
    # contrasted with blue or purplish-blue legs and feet.
    correct_answer_key = 'E'

    explanation = f"""
The Oophaga pumilio morph endemic to Isla Colón in the Bocas del Toro Archipelago is one of the most iconic. Its phenotype is distinguished by a bright red or orange-red body, which serves as a potent aposematic (warning) signal to predators. This vibrant body coloration is sharply contrasted by its limbs. The legs and feet of this morph are a distinct blue or sometimes purplish-blue.

Comparing this information with the given options:
- A, C, D, F, G, H, I are incorrect as they describe different color patterns or other morphs entirely.
- B is close but 'splotches' is less accurate than the general coloration of the legs.
- E provides the most accurate and widely recognized description.

Therefore, the best description for the phenotype of the 'Isla Colón' morph is:
    """

    print(textwrap.dedent(explanation))
    print(f"Correct Answer Choice ({correct_answer_key}): {answer_choices[correct_answer_key]}")

find_phenotype_description()