def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the correct physical exam maneuver.
    The code will print the step-by-step reasoning.
    """

    # --- Step 1: Analyze Patient Presentation ---
    print("Step 1: Analyzing the patient's key symptoms and history.")
    symptom = "Right lower extremity pain in an L4-S1 distribution."
    aggravating_factor = "Pain is intensified by lying supine."
    red_flag_history_1 = "Patient is on prednisone, causing immunosuppression."
    red_flag_history_2 = "History of ovarian cancer."
    print(f"- Main Symptom: {symptom}")
    print(f"- Key Feature: {aggravating_factor}. This is a 'red flag' that points away from simple mechanical back pain.")
    print(f"- Red Flags from History: {red_flag_history_1} and {red_flag_history_2}.")
    print("-" * 20)

    # --- Step 2: Formulate a Differential Diagnosis ---
    print("Step 2: Formulating a differential diagnosis based on red flags.")
    print("- The combination of pain worse when lying down, immunosuppression, and a history of cancer strongly suggests an inflammatory, infectious, or neoplastic process in the retroperitoneal space.")
    print("- A common diagnosis that fits this picture is a Psoas Abscess or a retroperitoneal tumor irritating the psoas muscle and the adjacent lumbosacral plexus nerve roots (L4-S1).")
    print("-" * 20)

    # --- Step 3: Identify Relevant Anatomy and Test ---
    print("Step 3: Identifying the relevant anatomy and the corresponding physical exam test.")
    print("- The psoas major muscle originates from the lumbar vertebrae and passes down through the pelvis.")
    print("- Irritation of this muscle is tested using the 'Psoas Sign'.")
    print("- The Psoas Sign is positive if pain is elicited by either passively stretching the psoas (hip extension) or actively contracting it against resistance (hip flexion).")
    print("-" * 20)
    
    # --- Step 4: Evaluate the Maneuver in the Question ---
    print("Step 4: Evaluating the specific maneuver described.")
    patient_position = "Left decubitus (lying on the left side), with the painful right leg on top."
    action = "The physician applies resistance to the extended right leg while the patient performs an action."
    print(f"- Patient Position: {patient_position}")
    print(f"- Action: {action}")
    print("- To test for a psoas process in this position, the patient would be asked to contract the psoas muscle by flexing the hip (bringing the knee to the chest) while the physician provides resistance.")
    print("-" * 20)

    # --- Step 5: Match the Test to the Answer Choices ---
    print("Step 5: Matching the correct action to the provided choices.")
    choices = {
        'A': 'Abduction',
        'B': 'Adduction',
        'C': 'Internal Rotation',
        'D': 'External Rotation',
        'E': 'Flexion',
        'F': 'Extension'
    }
    correct_action = "Flexion of the hip against resistance."
    correct_choice_letter = 'E'
    print(f"- The action that contracts the psoas muscle to elicit the psoas sign is hip flexion.")
    print(f"- Therefore, asking the patient to perform '{choices[correct_choice_letter]}' of the hip against resistance will confirm the diagnosis.")
    print("\nFinal Answer Calculation:")
    print(f"The correct maneuver is {choices[correct_choice_letter]}, which corresponds to choice {correct_choice_letter}.")


solve_clinical_case()