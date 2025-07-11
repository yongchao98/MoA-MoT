import sys

def solve_clinical_case():
    """
    This function analyzes the clinical vignette and determines the correct physical exam maneuver.
    """
    
    patient_presentation = {
        "Age": 55,
        "Sex": "Female",
        "Complaint": "Pain in lower right extremity, L4-S1 distribution",
        "Position for test": "Left decubitus (lying on left side)",
        "Limb being tested": "Right leg (extended)"
    }

    # The patient's symptoms (L4-S1 radiculopathy) suggest either a lumbar spine issue or
    # a nerve compression issue like piriformis syndrome.
    # The physical exam described is designed to differentiate these conditions.
    
    # The piriformis muscle's primary action on an EXTENDED hip is EXTERNAL ROTATION.
    # A provocative test for piriformis syndrome involves resisting this action
    # to see if contracting the muscle reproduces the patient's sciatic pain.

    # Let's analyze the options:
    # A. Abduction: Engages piriformis, but not its primary function in this position.
    # B. Adduction: Tests inner thigh muscles. Irrelevant.
    # C. Internal Rotation: Stretches the piriformis; the prompt mentions applying resistance to an action.
    # D. External Rotation: Directly tests the primary function of the piriformis in an extended hip. Contraction can compress the sciatic nerve.
    # E. Flexion: Tests hip flexors.
    # F. Extension: Tests gluteus maximus.

    conclusion = "Applying resistance to external rotation is the most specific maneuver to elicit pain from an irritated piriformis muscle, thereby confirming piriformis syndrome as the cause of the patient's sciatica-like symptoms."
    
    correct_answer_choice = "D"

    print("Clinical Reasoning:")
    print("-----------------")
    print(f"The patient's symptoms are consistent with sciatica (L4-S1 radiculopathy).")
    print("A key differential diagnosis for sciatica, especially with unremarkable imaging, is piriformis syndrome.")
    print("The physician is performing a provocative test in the left decubitus position on the extended right leg.")
    print("The primary function of the piriformis muscle on an extended hip is external rotation.")
    print("Therefore, to test for piriformis syndrome, the physician would resist this primary action.")
    print("\nConclusion:")
    print(conclusion)
    print(f"\nThe correct answer is {correct_answer_choice}: External Rotation.")

if __name__ == "__main__":
    # In a real scenario, we might process some data, but here we just explain the logic.
    # The prompt asked for code, so we use Python's print function to provide the step-by-step analysis and final answer.
    solve_clinical_case()