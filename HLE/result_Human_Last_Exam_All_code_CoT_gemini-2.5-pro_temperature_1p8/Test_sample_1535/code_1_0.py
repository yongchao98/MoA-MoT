def diagnose_rash_location():
    """
    Analyzes patient symptoms to determine the most likely location of a related rash.
    """
    # Key symptoms and findings from the patient's case
    patient_findings = {
        "muscle_weakness": True,
        "periorbital_erythema": True,  # A rash is already present on the eyelids
    }

    # Answer choices mapping to anatomical regions
    answer_choices = {
        "A": "Dorsum of the hands",
        "B": "Nose",
        "C": "Eyelids",
        "D": "Groin",
        "E": "Shoulders",
    }

    print("Step 1: Analyzing the patient's key symptoms.")
    print(f"- Proximal Muscle Weakness: Present")
    print(f"- Periorbital Erythema (rash on the eyelids): Present")
    print("\nStep 2: Forming a diagnosis.")
    print("The combination of muscle weakness and a characteristic rash around the eyes (Heliotrope rash) strongly suggests Dermatomyositis.")

    print("\nStep 3: Identifying other characteristic signs of Dermatomyositis.")
    print("Dermatomyositis is associated with several rashes. The most specific ones are:")
    print("- Heliotrope rash (on the Eyelids)")
    print("- Gottron's sign/papules (on the Dorsum of the hands)")
    print("- Shawl sign (on the Shoulders and back)")

    print("\nStep 4: Determining the expected location for another rash.")
    print("The patient already presents with a rash on the Eyelids (Periorbital Erythema).")
    print("A physician would next look for the other highly specific signs of Dermatomyositis.")
    print("Gottron's sign on the dorsum of the hands is one of the most specific findings.")
    
    final_answer_key = "A"
    final_answer_region = answer_choices[final_answer_key]
    
    print("\n--- Final Conclusion ---")
    print(f"Logical Equation: Dermatomyositis diagnosis + existing Eyelid rash -> strong expectation of rash on {final_answer_region}")
    print(f"The correct choice is '{final_answer_key}', which corresponds to the '{final_answer_region}'.")

diagnose_rash_location()
<<<A>>>