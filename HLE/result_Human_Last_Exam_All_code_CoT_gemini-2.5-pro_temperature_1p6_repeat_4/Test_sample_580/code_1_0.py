def solve_clinical_case():
    """
    This function analyzes the clinical scenario to identify the correct diagnostic maneuver.
    The reasoning is based on medical knowledge of anatomy and physical examination techniques.
    """
    
    # Step 1: Identify the primary symptom complex.
    symptoms = "Pain in the lower right extremity L4-S1 distribution."
    likely_condition = "Sciatica"
    
    # Step 2: Consider a key differential diagnosis that fits the context.
    differential_diagnosis = "Piriformis Syndrome"
    
    # Step 3: Recall the anatomy and function of the muscle involved.
    muscle = "Piriformis"
    primary_action = "External Rotation of the hip"
    
    # Step 4: Determine the principle of the provocative test.
    # To test for piriformis syndrome, the muscle is stretched to see if it compresses the sciatic nerve.
    # Stretching a muscle involves performing the action opposite to its primary function.
    provocative_maneuver = "The opposite of External Rotation"
    
    # Step 5: Identify the correct action from the given choices.
    # The opposite of External Rotation is Internal Rotation.
    correct_action = "Internal Rotation"
    answer_choice = "C"

    print("Clinical Reasoning:")
    print(f"1. The patient's symptoms ({symptoms}) point towards a diagnosis of {likely_condition}.")
    print(f"2. A key differential diagnosis is {differential_diagnosis}, where the sciatic nerve is compressed by the {muscle} muscle.")
    print(f"3. The primary function of the {muscle} muscle is {primary_action}.")
    print(f"4. A provocative test for this condition involves stretching the muscle by performing the opposite action.")
    print(f"5. Therefore, the maneuver that will stretch the piriformis and confirm the diagnosis is {correct_action}.")
    print("-" * 20)
    print(f"The correct answer choice is: {answer_choice}")

solve_clinical_case()