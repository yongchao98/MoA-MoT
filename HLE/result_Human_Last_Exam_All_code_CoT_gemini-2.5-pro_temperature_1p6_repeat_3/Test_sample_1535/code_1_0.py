def solve_medical_case():
    """
    Analyzes the clinical vignette to determine the correct anatomical region for the rash.
    """
    # Step 1: Define key signs from the vignette
    muscle_symptoms = ["muscle weakness", "myalgia"]
    skin_symptom = "periorbital erythema" # This means redness around the eyes.

    # Step 2: Combine symptoms to reach a likely diagnosis
    # Muscle inflammation (Myositis) + Skin inflammation (Dermatitis) points to Dermatomyositis.
    diagnosis = "Dermatomyositis"

    # Step 3: Relate the specific sign to an anatomical location
    # In Dermatomyositis, "periorbital erythema" is the clinical term for a
    # classic rash called the Heliotrope rash.
    rash_name = "Heliotrope Rash"
    rash_location = "Eyelids"

    # Step 4: Print the logical deduction (the "equation")
    print("Logical Analysis:")
    print(f"1. The patient's muscle symptoms ({', '.join(muscle_symptoms)}) suggest Myositis.")
    print(f"2. The patient's skin symptom ('{skin_symptom}') suggests Dermatitis.")
    print(f"3. Myositis + Dermatitis points to the diagnosis of {diagnosis}.")
    print(f"4. The term '{skin_symptom}' describes the {rash_name}, which is characteristically located on the {rash_location}.")
    
    print("\nFinal Logical Equation:")
    print(f"'{skin_symptom}' == 'Rash on the {rash_location}'")
    print(f"Therefore, the correct anatomical region is {rash_location}.")
    print("Answer choice C corresponds to this location.")

solve_medical_case()