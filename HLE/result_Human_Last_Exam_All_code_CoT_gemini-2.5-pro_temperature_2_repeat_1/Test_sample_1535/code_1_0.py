def solve_clinical_case():
    """
    Analyzes clinical findings to determine the expected location of a rash.
    """
    # Step 1: Define the key clinical findings from the patient's case.
    patient_findings = {
        "chief_complaint": "Post-operational back pain",
        "history": ["muscle weakness", "myalgia", "arthralgia", "congestive heart disease"],
        "physical_exam": ["pale face", "periorbital recession", "periorbital erythema"]
    }

    # Step 2: Identify the most significant signs for diagnosis.
    key_symptom_muscle = "muscle weakness"
    key_symptom_skin = "periorbital erythema" # Redness around the eyes.

    print("--- Clinical Reasoning ---")
    print(f"1. The patient presents with a significant finding of '{key_symptom_muscle}'.")
    print(f"2. The physical exam reveals '{key_symptom_skin}'.")
    print("")

    # Step 3 & 4: Link the signs to a specific diagnosis and its characteristic rash.
    diagnosis = "Dermatomyositis"
    characteristic_rash_name = "Heliotrope rash"
    rash_location = "Eyelids"
    
    print("--- Diagnosis ---")
    print(f"The combination of muscle weakness and a rash around the eyes is a classic presentation of {diagnosis}.")
    print(f"The sign '{key_symptom_skin}' describes the {characteristic_rash_name}, a hallmark of this condition.")
    print("")
    
    # Step 5: Determine the anatomical location and select the correct answer.
    print("--- Conclusion ---")
    print(f"The final equation for the answer is: '{key_symptom_skin}' points to '{characteristic_rash_name}' which is located on the '{rash_location}'.")
    print(f"Based on this, the rash is expected on the anatomical region listed in choice C.")

solve_clinical_case()