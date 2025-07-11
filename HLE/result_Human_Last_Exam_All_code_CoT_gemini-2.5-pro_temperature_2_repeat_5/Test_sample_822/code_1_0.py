def solve_clinical_case():
    """
    Analyzes the clinical findings to determine the most likely diagnosis.
    This mimics the process of adding up clinical clues to reach a conclusion.
    """
    # Key clinical findings from the case vignette are treated as the "numbers".
    age = 62
    smoking_history_pack_years = 20
    symptom_1 = "Polyarthritis (pain in wrists, ankles, elbows)"
    symptom_2 = "Multiple pulmonary nodules"
    symptom_3 = "Systemic symptoms (fatigue, confusion, bruising)"
    symptom_4 = "Cutaneous lesions"
    symptom_5 = "Symptoms suggestive of upper airway/sinus involvement (difficulty swallowing)"

    print("Analyzing the key clinical findings ('the numbers'):")
    print(f"1. Age: {age}")
    print(f"2. Smoking History (Pack-Years): {smoking_history_pack_years}")
    print(f"3. Musculoskeletal Finding: {symptom_1}")
    print(f"4. Primary Radiologic Finding: {symptom_2}")
    print(f"5. Other Systemic Findings: {symptom_3}")
    print(f"6. Dermatologic Finding: {symptom_4}")
    print(f"7. Head and Neck Finding: {symptom_5}")
    
    print("\nFormulating the Diagnostic 'Equation':")
    print(f"'{symptom_1}' + '{symptom_2}' + '{symptom_5}'")
    print("This combination strongly suggests a disease process affecting the joints, lungs, and upper airways.")
    
    print("\nAdding further evidence:")
    print("When combined with systemic, neurological, and cutaneous symptoms, the diagnosis points towards a multi-system small-vessel vasculitis.")
    
    diagnosis = "Granulomatosis with Polyangiitis (GPA)"
    
    print("\nFinal Conclusion:")
    print(f"The constellation of findings is most characteristic of {diagnosis}.")
    print("The patient's immunocompromised state (from steroid treatment for the disease) led to a fatal superimposed infection.")

solve_clinical_case()