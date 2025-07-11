import collections

def solve_clinical_vignette():
    """
    This script analyzes the clinical findings from the vignette to determine the most likely location of a rash.
    """
    # Key findings from the patient's physical exam and history
    patient_findings = {
        "muscle weakness",
        "myalgia",
        "periorbital erythema" # Redness around the eyes
    }

    # Database of diseases and their pathognomonic signs
    # Using a set for signs allows for easy checking of subsets
    disease_database = {
        "Dermatomyositis": {"muscle weakness", "periorbital erythema"},
        "Systemic Lupus Erythematosus": {"malar rash", "arthralgia"},
        "Psoriatic Arthritis": {"psoriasis", "arthralgia"}
    }

    # Database of diseases and their associated rash locations
    rash_locations = {
        "Dermatomyositis": [
            ("Heliotrope rash", "Eyelids"),
            ("Gottron's papules", "Dorsum of the hands"),
            ("Shawl sign", "Shoulders")
        ]
    }

    # Step 1: Determine the most likely diagnosis
    likely_diagnosis = None
    for disease, key_signs in disease_database.items():
        if key_signs.issubset(patient_findings):
            likely_diagnosis = disease
            break

    # Step 2: Explain the reasoning and find the answer
    print("Clinical Reasoning Steps:")
    print("-------------------------")
    if likely_diagnosis:
        print(f"1. The combination of 'muscle weakness' and 'periorbital erythema' strongly suggests a diagnosis of {likely_diagnosis}.")

        # Find the locations associated with the diagnosis
        locations = rash_locations.get(likely_diagnosis, [])
        if locations:
            print(f"2. {likely_diagnosis} is associated with several characteristic rashes:")
            for rash_name, location in locations:
                print(f"   - {rash_name} (on the {location})")

            print("\n3. The specific physical exam finding 'periorbital erythema' (redness around the eyes) directly corresponds to the Heliotrope rash.")
            print("   The anatomical location for the Heliotrope rash is the Eyelids.")

            # Match the location to the answer choices
            answer_choices = {
                'A': "Dorsum of the hands",
                'B': "Nose",
                'C': "Eyelids",
                'D': "Groin",
                'E': "Shoulders"
            }
            final_answer_letter = None
            for letter, description in answer_choices.items():
                if description == "Eyelids":
                    final_answer_letter = letter
                    break
            
            print("\nConclusion:")
            print(f"Based on this reasoning, the expected anatomical region for the rash is '{answer_choices[final_answer_letter]}'.")
            print(f"The correct option is: {final_answer_letter}")
        else:
            print("No characteristic rash locations found for the diagnosis in the database.")
    else:
        print("Could not determine a likely diagnosis based on the provided findings.")

# Run the analysis
solve_clinical_vignette()