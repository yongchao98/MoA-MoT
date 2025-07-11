def solve_clinical_case():
    """
    This script analyzes a clinical case to determine the expected location of a rash
    by modeling the diagnostic process.
    """
    # Step 1: Define the patient's key clinical findings from the case vignette.
    # The combination of muscle weakness and periorbital erythema is highly suggestive.
    print("Step 1: Parsing the patient's key clinical findings.")
    patient_findings = {
        "muscle weakness",
        "myalgia",             # muscle pain
        "arthralgia",          # joint pain
        "periorbital erythema" # redness around the eyes, i.e., Heliotrope rash
    }
    print(f"-> Key findings identified: muscle weakness, myalgia, and periorbital erythema.\n")

    # Step 2: Define profiles for the most relevant differential diagnoses.
    # Periorbital erythema (Heliotrope rash) + muscle weakness is classic for Dermatomyositis.
    print("Step 2: Comparing findings against disease profiles.")
    disease_profiles = {
        "Dermatomyositis": {
            "symptoms": {"muscle weakness", "myalgia", "arthralgia"},
            "rashes": {"heliotrope rash", "gottron's sign"},
            "rash_locations": {
                "heliotrope rash": "Eyelids",
                "gottron's sign": "Dorsum of the hands",
                "shawl sign": "Shoulders"
            }
        },
        "Systemic Lupus Erythematosus": {
            "symptoms": {"arthralgia"},
            "rashes": {"malar rash"},
            "rash_locations": {
                "malar rash": "Nose"
            }
        }
    }
    
    # Step 3: Identify the most likely diagnosis.
    # In this case, the combination of a specific rash and muscle symptoms points directly to one diagnosis.
    patient_has_heliotrope = "periorbital erythema" in patient_findings
    patient_has_myopathy = "muscle weakness" in patient_findings

    likely_diagnosis = ""
    if patient_has_heliotrope and patient_has_myopathy:
        likely_diagnosis = "Dermatomyositis"
        
    print(f"-> The combination of myopathy and a heliotrope rash strongly indicates a diagnosis of {likely_diagnosis}.\n")

    # Step 4: Identify other expected skin findings for the diagnosed condition.
    # The question asks for an *expected* rash, implying another characteristic site.
    print(f"Step 4: Identifying other characteristic rashes for {likely_diagnosis}.")
    
    # The other pathognomonic (highly specific) rash for Dermatomyositis is Gottron's sign.
    expected_rash_name = "gottron's sign"
    expected_location = disease_profiles[likely_diagnosis]["rash_locations"][expected_rash_name]

    print(f"-> The patient already presents with a Heliotrope rash on the Eyelids.")
    print(f"-> The other classic, expected rash is {expected_rash_name}, which appears on the '{expected_location}'.\n")

    # Step 5: Match the anatomical location to the given answer choices.
    print("Step 5: Matching the location to the answer choices.")
    answer_choices = {
        'A': "Dorsum of the hands",
        'B': "Nose",
        'C': "Eyelids",
        'D': "Groin",
        'E': "Shoulders"
    }

    final_answer_letter = ""
    for choice, location in answer_choices.items():
        if location == expected_location:
            final_answer_letter = choice
            break
            
    print(f"-> The location '{expected_location}' corresponds to answer choice: {final_answer_letter}.")

solve_clinical_case()