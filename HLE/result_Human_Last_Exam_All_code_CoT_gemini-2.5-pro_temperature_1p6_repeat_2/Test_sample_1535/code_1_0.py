def solve_medical_case():
    """
    This function analyzes the clinical vignette to determine the most likely location for a rash.
    """
    # 1. Define the patient's key findings from the vignette.
    patient_findings = {
        "muscle weakness",
        "myalgia",
        "arthralgia",
        "periorbital erythema"  # This is the most specific physical finding.
    }

    # 2. Define the classic features of the suspected diagnosis.
    # The combination of muscle inflammation and skin rashes points to Dermatomyositis.
    dermatomyositis_rash_locations = {
        "Heliotrope rash": "Eyelids",
        "Gottron's sign/papules": "Dorsum of the hands",
        "Shawl sign": "Shoulders",
        "V sign": "Neck and upper chest"
    }

    # 3. Explain the reasoning process.
    print("Step 1: Analyzing the patient's presentation.")
    print("The patient exhibits muscle weakness, myalgia (muscle pain), and a key physical finding of 'periorbital erythema' (redness around the eyes).")
    print("\nStep 2: Identifying the underlying condition.")
    print("This combination of symptoms strongly suggests the diagnosis of Dermatomyositis, an inflammatory disease affecting the muscles and skin.")
    print("\nStep 3: Correlating the diagnosis with the physical finding.")
    print("In Dermatomyositis, 'periorbital erythema' is the clinical term for a characteristic rash called the 'Heliotrope rash'.")
    
    # 4. Determine the anatomical location.
    heliotrope_rash_location = dermatomyositis_rash_locations["Heliotrope rash"]
    print(f"The Heliotrope rash is located specifically on the '{heliotrope_rash_location}'.")

    # 5. Evaluate the answer choices.
    answer_choices = {
        'A': 'Dorsum of the hands',
        'B': 'Nose',
        'C': 'Eyelids',
        'D': 'Groin',
        'E': 'Shoulders'
    }
    
    correct_answer_key = 'C'
    
    print("\nStep 4: Conclusion.")
    print(f"The question asks which region is expected to have a rash. Since the patient has 'periorbital erythema', a rash is, by definition, present on the {answer_choices[correct_answer_key]}.")
    print(f"Therefore, the correct answer is C.")

solve_medical_case()