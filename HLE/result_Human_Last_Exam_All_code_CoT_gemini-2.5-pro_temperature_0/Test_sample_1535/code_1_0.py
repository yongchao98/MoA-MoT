def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the location of a rash.
    """
    # Step 1: Identify the key clinical features from the patient's presentation.
    # The combination of systemic symptoms like muscle weakness, myalgia, and fatigue,
    # along with a specific skin finding, is crucial.
    print("Step 1: Analyzing the patient's key signs and symptoms.")
    symptoms = ["muscle weakness", "myalgia", "fatigue"]
    skin_finding = "periorbital erythema"
    print(f"The patient has systemic symptoms ({', '.join(symptoms)}) and a key skin finding: '{skin_finding}'.")
    print("This constellation of findings is highly suggestive of Dermatomyositis.")
    print("-" * 40)

    # Step 2: Interpret the specific physical exam finding.
    # "Periorbital" refers to the area around the eyes. In Dermatomyositis, this
    # rash is known as the "Heliotrope rash".
    print("Step 2: Interpreting the physical exam finding.")
    print(f"The term '{skin_finding}' describes redness around the eye sockets.")
    print("This is the clinical description for the 'Heliotrope rash', a pathognomonic (highly specific) sign of Dermatomyositis.")
    print("-" * 40)

    # Step 3: Determine the precise anatomical location of this rash.
    # The Heliotrope rash is characteristically found on the eyelids.
    location = "Eyelids"
    print("Step 3: Identifying the anatomical location.")
    print(f"The Heliotrope rash is located on the upper {location}.")
    print("-" * 40)

    # Step 4: Match the location to the given answer choices.
    answer_choices = {
        'A': 'Dorsum of the hands',
        'B': 'Nose',
        'C': 'Eyelids',
        'D': 'Groin',
        'E': 'Shoulders'
    }
    correct_choice_letter = 'C'
    print("Step 4: Selecting the correct answer.")
    print(f"The clinical finding of 'periorbital erythema' directly corresponds to a rash on the '{answer_choices[correct_choice_letter]}'.")
    print(f"Therefore, the correct answer is C.")

# Execute the reasoning process
solve_clinical_case()
<<<C>>>