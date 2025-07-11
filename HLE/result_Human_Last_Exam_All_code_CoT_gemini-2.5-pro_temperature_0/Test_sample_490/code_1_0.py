def find_phenotype_description():
    """
    This function identifies and explains the correct phenotype for the
    "Isla Colón" morph of Oophaga pumilio.
    """
    question = 'Which of the following best describes the phenotype of the "Isla Colón" morph of Oophaga pumilio?'
    
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

    correct_answer_key = 'F'
    
    explanation = (
        "The 'Isla Colón' morph of the Oophaga pumilio poison dart frog, "
        "found in the Bocas del Toro Archipelago of Panamá, is famously characterized "
        "by its bright, solid red or orange-red coloration. While other morphs, "
        "like the well-known 'blue-jeans' frog (often from other locations), "
        "have blue legs, the Isla Colón variety is typically monochromatic red."
    )

    print(f"Question: {question}\n")
    print("Answer Choices:")
    for key, value in answer_choices.items():
        print(f"{key}. {value}")
    
    print("\n----------------------------------------\n")
    print(f"Correct Answer: {correct_answer_key}")
    print(f"Description: {answer_choices[correct_answer_key]}")
    print(f"\nExplanation: {explanation}")

find_phenotype_description()