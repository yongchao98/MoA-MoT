def solve_clinical_vignette():
    """
    This function analyzes a clinical vignette to determine the expected location of a rash.
    """
    # Step 1: Define the key clinical findings from the patient's presentation.
    patient_findings = {
        "Primary Complaint": "Post-operational back pain",
        "Systemic Symptoms": ["Fatigue", "Arthralgia", "Myalgia", "Muscle weakness"],
        "Physical Exam": ["Periorbital recession", "Periorbital erythema"]
    }

    # Step 2: Identify the most indicative sign and the likely diagnosis.
    key_finding = "Periorbital erythema"
    likely_diagnosis = "Dermatomyositis"
    
    print("Clinical Reasoning Steps:")
    print("-" * 25)
    
    # Step 3: Connect the key finding to its specific medical term and location.
    print(f"1. The patient presents with muscle weakness and a key physical exam finding: '{key_finding}'.")
    
    characteristic_rash_name = "Heliotrope rash"
    print(f"2. The combination of muscle weakness and '{key_finding}' is a classic presentation for {likely_diagnosis}.")
    
    rash_location = "Eyelids"
    print(f"3. In {likely_diagnosis}, '{key_finding}' is the primary component of the pathognomonic (highly specific) '{characteristic_rash_name}'.")

    # Step 4: Determine the answer based on the location of the characteristic rash.
    print(f"4. The '{characteristic_rash_name}' is located on the {rash_location}.")
    
    answer_choices = {
        'A': 'Dorsum of the hands',
        'B': 'Nose',
        'C': 'Eyelids',
        'D': 'Groin',
        'E': 'Shoulders'
    }

    final_answer_code = 'C'
    
    print("-" * 25)
    print(f"Conclusion: The rash is expected on the {answer_choices[final_answer_code]}.")
    print(f"This corresponds to answer choice {final_answer_code}.")

solve_clinical_vignette()