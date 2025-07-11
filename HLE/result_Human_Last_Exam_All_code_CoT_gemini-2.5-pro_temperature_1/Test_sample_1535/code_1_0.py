def solve_medical_case():
    """
    This function analyzes the clinical vignette to determine the expected location of a rash.
    """

    # Step 1: Isolate the key clinical findings from the case description.
    key_finding = "Periorbital erythema (redness around the eyes)"
    associated_symptoms = ["Muscle weakness", "Myalgia", "Arthralgia"]

    print("Clinical Reasoning Steps:")
    print("========================")

    # Step 2: Define the medical term for the key finding and its significance.
    medical_term = "Heliotrope rash"
    associated_disease = "Dermatomyositis"
    
    print(f"1. The patient's most specific physical finding is: '{key_finding}'.")
    print(f"   + This is a classic sign known as a '{medical_term}'.")

    # Step 3: Connect the sign to the most likely diagnosis.
    print(f"2. A '{medical_term}' combined with symptoms like '{associated_symptoms[0]}' and '{associated_symptoms[1]}' is highly characteristic of the condition '{associated_disease}'.")

    # Step 4: Determine the anatomical location based on the definition of the sign.
    location_of_heliotrope_rash = "Eyelids"
    answer_choice_letter = "C"
    
    print(f"3. By definition, a '{medical_term}' is a rash located on the '{location_of_heliotrope_rash}'.")

    # Final Conclusion
    print("\n--- Conclusion ---")
    print(f"The question asks for the expected location of the rash. Based on the finding of 'periorbital erythema', the answer is:")
    print(f"Final Answer = {answer_choice_letter}. {location_of_heliotrope_rash}")

solve_medical_case()