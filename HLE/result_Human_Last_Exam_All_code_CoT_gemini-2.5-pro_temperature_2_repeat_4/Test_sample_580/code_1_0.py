import sys
import io

# Redirect stdout to a string buffer to capture all printed output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_clinical_vignette():
    """
    This script analyzes a clinical vignette to determine the correct diagnostic maneuver.
    """
    
    # Step 1: Define the key clinical findings from the vignette.
    patient_profile = {
        "Complaint": "Pain in the lower right extremity, following the L4-S1 dermatome distribution.",
        "Key Finding 1": "Unremarkable X-ray imaging, suggesting a non-bony or soft tissue etiology.",
        "Key Finding 2": "Pain is intensified by lying supine, which is atypical for classic lumbar radiculopathy.",
        "Examination Position": "Patient is in the left decubitus position (lying on the left side).",
        "Test Setup": "Physician applies resistance to the extended right leg."
    }

    # Step 2: Formulate the most likely diagnosis based on the evidence.
    likely_diagnosis = "Piriformis Syndrome"
    reasoning_for_diagnosis = "Piriformis syndrome, where the piriformis muscle compresses the sciatic nerve, is a primary consideration as it mimics sciatica but originates from a soft tissue structure in the buttock."

    print("--- Diagnostic Analysis ---")
    print(f"1. The patient's symptoms ({patient_profile['Complaint']}) are consistent with sciatica.")
    print(f"2. {patient_profile['Key Finding 1']} and {patient_profile['Key Finding 2']}")
    print(f"3. Based on these findings, the most probable diagnosis is: {likely_diagnosis}.")
    
    # Step 3: Determine the confirmatory physical exam maneuver.
    piriformis_function = "The primary function of the piriformis muscle is EXTERNAL ROTATION of the hip."
    test_logic = "To confirm piriformis syndrome, a provocative test involves having the patient contract the muscle against resistance to see if it reproduces the sciatic nerve pain."
    
    answer_choices = {
        "A": "Abduction",
        "B": "Adduction",
        "C": "Internal Rotation",
        "D": "External Rotation",
        "E": "Flexion",
        "F": "Extension"
    }
    
    correct_action = "External Rotation"
    correct_choice_letter = "D"
    
    print("\n--- Confirmatory Test Rationale ---")
    print(f"Function: {piriformis_function}")
    print(f"Logic: {test_logic}")
    print(f"Conclusion: In the described position ({patient_profile['Examination Position']}), applying resistance to {correct_action} of the right leg is the most specific test to confirm {likely_diagnosis}.")
    
    # Step 4: Output the final answer from the choices provided.
    print(f"\nThe correct action is '{answer_choices[correct_choice_letter]}'.")


# Execute the diagnostic function
solve_clinical_vignette()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the console
print(output)

# Print the final answer in the required format
print("<<<D>>>")