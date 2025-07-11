def solve_clinical_case():
    """
    This script analyzes a clinical vignette to determine the expected location of a rash.
    """
    # Step 1: Define the key clinical findings from the case.
    key_symptoms = ["muscle weakness", "myalgia", "arthralgia", "fatigue"]
    key_physical_finding = "periorbital erythema (redness around the eyes)"

    print("Step-by-step Clinical Reasoning:\n")
    print(f"1. The patient presents with a combination of muscle-related issues ({', '.join(key_symptoms)}) and a key skin finding: {key_physical_finding}.")

    # Step 2: Correlate findings with a likely diagnosis.
    diagnosis = "Dermatomyositis"
    print(f"\n2. This constellation of symptoms is highly characteristic of an inflammatory myopathy known as {diagnosis}.")

    # Step 3: Identify the classic skin manifestation associated with the diagnosis and findings.
    classic_rash_name = "Heliotrope rash"
    classic_rash_location = "Eyelids"
    print(f"\n3. The finding of '{key_physical_finding}' is the classic description of a '{classic_rash_name}'.")

    # Step 4: Determine the specific anatomical location of this rash.
    print(f"\n4. The {classic_rash_name} is specifically located on the '{classic_rash_location}'.")
    
    # Step 5: Final conclusion based on the provided answer choices.
    answer_choices = {
        'A': 'Dorsum of the hands',
        'B': 'Nose',
        'C': 'Eyelids',
        'D': 'Groin',
        'E': 'Shoulders'
    }
    
    final_answer_letter = 'C'
    final_answer_description = answer_choices[final_answer_letter]

    print(f"\nConclusion: Based on the clinical picture pointing to {diagnosis}, the 'periorbital erythema' corresponds to the Heliotrope rash. Therefore, the anatomical region expected to have the rash is the {final_answer_description}.")


solve_clinical_case()
<<<C>>>