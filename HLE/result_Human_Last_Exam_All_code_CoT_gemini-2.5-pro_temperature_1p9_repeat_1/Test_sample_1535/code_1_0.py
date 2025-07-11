def analyze_medical_case():
    """
    This script analyzes a clinical vignette to determine the location of a characteristic rash.
    """
    # Define the key symptoms and signs from the patient's case.
    symptoms = {
        "Musculoskeletal": ["myalgia", "muscle weakness", "arthralgia"],
        "Systemic": ["fatigue"],
        "Dermatologic": ["periorbital recession and erythema"]
    }

    # Step 1: Synthesize findings to form a likely diagnosis.
    # The combination of muscle inflammation (myositis) and skin inflammation (dermatitis)
    # is the hallmark of Dermatomyositis.
    diagnosis = "Dermatomyositis"
    print(f"Analysis Step 1: The patient's constellation of symptoms (muscle weakness, fatigue, skin findings) strongly suggests a diagnosis of {diagnosis}.")

    # Step 2: Interpret the specific physical exam finding.
    # "Periorbital erythema" means redness around the eye sockets.
    # In Dermatomyositis, this is a classic finding known as the "heliotrope rash".
    finding = "periorbital erythema"
    interpretation = "heliotrope rash"
    print(f"Analysis Step 2: The finding of '{finding}' corresponds to the pathognomonic '{interpretation}'.")

    # Step 3: Identify the anatomical location of the key sign.
    # The heliotrope rash is characteristically located on the eyelids.
    location = "Eyelids"
    print(f"Analysis Step 3: The {interpretation} is located on the {location}.")

    # Step 4: Evaluate the provided answer choices.
    answer_choices = {
        "A": "Dorsum of the hands",
        "B": "Nose",
        "C": "Eyelids",
        "D": "Groin",
        "E": "Shoulders"
    }
    
    # Conclusion
    # While rashes on the dorsum of the hands (A) and shoulders (E) also occur in Dermatomyositis,
    # the clinical description explicitly points to the periorbital area.
    correct_choice = "C"
    
    print("\nConclusion:")
    print(f"The clinical sign described ('{finding}') directly points to the anatomical region in choice {correct_choice}.")
    print(f"Therefore, the expected region for the rash is: {answer_choices[correct_choice]}.")

# Execute the analysis.
analyze_medical_case()