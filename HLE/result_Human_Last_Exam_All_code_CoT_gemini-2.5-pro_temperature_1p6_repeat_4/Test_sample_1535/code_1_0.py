import textwrap

def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the location of a characteristic rash.
    """
    # Key finding from the physical exam in the vignette
    key_finding = "periorbital erythema"

    # Medical knowledge base linking findings to anatomical locations
    medical_definitions = {
        "periorbital erythema": {
            "common_name": "Heliotrope rash",
            "location": "Eyelids",
            "associated_condition": "Dermatomyositis"
        }
    }

    # Answer choices provided
    answer_choices = {
        "A": "Dorsum of the hands",
        "B": "Nose",
        "C": "Eyelids",
        "D": "Groin",
        "E": "Shoulders"
    }

    # Step-by-step reasoning
    print("Step 1: Identify the most specific clinical sign from the physical exam.")
    print(f"The key finding is: '{key_finding}'.\n")

    print("Step 2: Define the clinical sign and its anatomical location.")
    if key_finding in medical_definitions:
        definition = medical_definitions[key_finding]
        explanation = (f"The term '{key_finding}' describes redness (erythema) in the area around the eyes (periorbital). "
                       f"This specific rash is also known as a Heliotrope rash and is located on the eyelids. It is a classic sign of {definition['associated_condition']}.")
        
        # Use textwrap for better formatting of long lines
        wrapper = textwrap.TextWrapper(width=80)
        print(wrapper.fill(text=explanation))
        print("\n")
        
        correct_location = definition['location']
        
        print(f"Step 3: Match the anatomical location with the given answer choices.")
        for choice, location in answer_choices.items():
            if location == correct_location:
                print(f"The location '{correct_location}' matches answer choice {choice}.")
                final_answer = choice
                break
    else:
        print("Key finding not in knowledge base.")
        final_answer = None

    print("\n--- Conclusion ---")
    if final_answer:
        print(f"The physical exam explicitly describes erythema around the eyes. This rash is located on the eyelids.")
        print(f"Therefore, the correct answer is {final_answer}.")
    else:
        print("Could not determine the answer.")

solve_clinical_case()