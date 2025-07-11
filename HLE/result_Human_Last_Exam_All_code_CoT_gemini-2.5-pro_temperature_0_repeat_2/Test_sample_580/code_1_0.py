def solve_clinical_case():
    """
    This function simulates the reasoning process to solve the clinical vignette.
    """
    # Step 1: Define the key information from the case.
    symptoms = "L4-S1 pain, unremarkable X-ray"
    exam_position = "Left decubitus (side-lying)"
    exam_action = "Resistance applied to the EXTENDED right leg"
    likely_diagnosis = "Piriformis Syndrome"

    print(f"Patient Presentation: {symptoms}")
    print(f"Physical Exam Setup: {exam_position}, {exam_action}")
    print(f"Based on the presentation (sciatic-like pain with normal imaging) and exam setup, the likely diagnosis is {likely_diagnosis}.")
    print("-" * 20)

    # Step 2: Analyze the anatomy and function relevant to the diagnosis.
    print("To confirm Piriformis Syndrome, we must test the piriformis muscle.")
    print("The function of the piriformis muscle changes with hip position:")
    print(" - When the hip is FLEXED, it acts as an abductor.")
    print(" - When the hip is EXTENDED, it acts as an EXTERNAL ROTATOR.")
    print("-" * 20)

    # Step 3: Determine the correct maneuver based on the exam setup.
    print("The case states the patient's leg is EXTENDED.")
    print("The test involves the patient performing an action against resistance.")
    print("Therefore, the correct maneuver to activate the piriformis muscle in this position is to have the patient perform its primary action.")
    
    final_action = "External Rotation"
    answer_choice = "D"

    print(f"The required action is: {final_action}")
    print("-" * 20)
    print(f"Final Answer Choice: {answer_choice}. {final_action}")


solve_clinical_case()