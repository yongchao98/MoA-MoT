import sys
import io

# Redirect stdout to capture the print statements
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_medical_case():
    """
    This function analyzes the patient's clinical case to identify the location of a characteristic rash.
    """
    # Step 1: Define the patient's key clinical findings from the case description.
    patient_findings = {
        "muscle weakness",
        "periorbital erythema"  # This is the key physical exam finding.
    }

    # Step 2: Define a simplified knowledge base of diseases and their key signs.
    # The combination of muscle weakness and a specific rash is crucial.
    disease_profiles = {
        "Dermatomyositis": {
            "core_symptoms": {"muscle weakness"},
            "rash_sign_map": {
                "periorbital erythema": "Eyelids",  # This is the Heliotrope rash
                "erythema over knuckles": "Dorsum of the hands" # This is Gottron's sign
            }
        },
        "Systemic Lupus Erythematosus": {
            "core_symptoms": {"arthralgia"},
            "rash_sign_map": {
                "malar erythema": "Nose"
            }
        }
    }

    # Step 3: Analyze the findings to determine the most likely diagnosis and rash location.
    print("Analyzing the clinical case:")
    print(f"Patient's key findings: {', '.join(patient_findings)}")

    likely_diagnosis = None
    rash_location = None

    # Find the disease that matches both the core symptom and the specific rash.
    for disease, profile in disease_profiles.items():
        has_core_symptom = not profile["core_symptoms"].isdisjoint(patient_findings)
        
        # Check if any of the patient's rash signs match the disease profile
        matching_rash_sign = None
        for sign in patient_findings:
            if sign in profile["rash_sign_map"]:
                matching_rash_sign = sign
                break

        if has_core_symptom and matching_rash_sign:
            likely_diagnosis = disease
            rash_location = profile["rash_sign_map"][matching_rash_sign]
            break

    # Step 4: Formulate the conclusion based on the analysis.
    if likely_diagnosis and rash_location:
        print(f"\nThe combination of '{list(patient_findings)[0]}' and '{list(patient_findings)[1]}' is characteristic of {likely_diagnosis}.")
        print(f"The specific finding of '{matching_rash_sign}' corresponds to a rash on the {rash_location}.")
    else:
        print("\nCould not determine a conclusive diagnosis from the provided key findings.")

    # Step 5: Match the location to the given answer choices.
    answer_choices = {
        "A": "Dorsum of the hands",
        "B": "Nose",
        "C": "Eyelids",
        "D": "Groin",
        "E": "Shoulders"
    }

    final_answer_letter = [letter for letter, desc in answer_choices.items() if desc == rash_location][0]

    print(f"\nTherefore, the expected anatomical region for the rash is the '{rash_location}'.")
    print(f"This corresponds to answer choice: {final_answer_letter}")
    return final_answer_letter

# Execute the analysis
final_answer = solve_medical_case()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())

# The final answer format required by the prompt
# print(f"<<<{final_answer}>>>")