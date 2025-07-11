def diagnose_skin_condition():
    """
    This script analyzes a clinical vignette to determine the most likely diagnosis.
    The reasoning is based on matching the patient's risk factors and physical exam
    findings with the characteristics of the differential diagnoses.
    """
    
    # Patient Data & Key Findings
    patient_age = 64
    patient_bmi = 39  # Obesity
    is_smoker = True
    
    # Lesion locations are key intertriginous areas
    lesion_locations = ["axillary folds", "inframammary folds", "inguinal regions"]
    
    # Lesion types include a mix, but one is highly specific
    lesion_types = {
        "axillary": "large bullae",
        "inframammary": "erythematous skin regions with plaques",
        "inguinal": "purulent nodules"
    }

    # Answer Choices
    diagnoses = {
        "A": "Malignant Intertrigo",
        "B": "Allergic contact dermatitis",
        "C": "Hidradenitis Supportiva",
        "D": "Atopic dermatitis",
        "E": "Psoriasis"
    }

    print("Clinical Reasoning:")
    print(f"1. Patient presents with a BMI of {patient_bmi} (obesity) and is a smoker. Both are major risk factors for Hidradenitis Suppurativa (HS).")
    print(f"2. The lesions are located in classic HS sites: {', '.join(lesion_locations)}.")
    print(f"3. The most specific finding is the presence of '{lesion_types['inguinal']}' in the groin.")
    print("4. This finding, in combination with the risk factors and locations, is hallmark for HS.")
    
    final_diagnosis_code = "C"
    final_diagnosis_name = diagnoses[final_diagnosis_code]

    print("\nConclusion:")
    print(f"The most likely diagnosis is {final_diagnosis_code}: {final_diagnosis_name}.")

diagnose_skin_condition()