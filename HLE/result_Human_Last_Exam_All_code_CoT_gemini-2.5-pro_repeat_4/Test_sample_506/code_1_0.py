def recommend_htn_medications():
    """
    Analyzes the patient case and recommends a 3-drug regimen for resistant hypertension.
    The code will print the recommendations and the clinical rationale for each choice.
    """

    # Medication recommendations are based on the patient's profile:
    # - Resistant Hypertension (Stage II)
    # - Comorbidities: Diabetes (A1C = 6.7%), African American ethnicity
    # - Specific medication intolerances noted in the case.

    medication_1 = {
        "name": "Chlorthalidone",
        "class": "Thiazide-like Diuretic",
        "rationale": "This is a first-line agent for hypertension, especially effective in African American patients and for resistant cases. It is a powerful, long-acting diuretic and is notably absent from the patient's list of contraindicated medications."
    }

    medication_2 = {
        "name": "Amlodipine",
        "class": "Dihydropyridine Calcium Channel Blocker",
        "rationale": "This is another first-line agent for hypertension in African American patients. It provides potent blood pressure lowering through vasodilation, offering a synergistic effect with a diuretic."
    }

    medication_3 = {
        "name": "Lisinopril",
        "class": "ACE Inhibitor",
        "rationale": "This is recommended due to the patient's compelling comorbidity of diabetes, indicated by her A1C of 6.7%. It provides essential kidney protection and further helps control blood pressure by targeting the RAAS pathway."
    }

    print("Recommendation for MM's Hypertension Treatment:")
    print("-" * 80)
    print("To maximize her hypertension treatment, the following 3-drug regimen is recommended:\n")

    print(f"1. Medication: {medication_1['name']} ({medication_1['class']})")
    print(f"   Rationale: {medication_1['rationale']}\n")

    print(f"2. Medication: {medication_2['name']} ({medication_2['class']})")
    print(f"   Rationale: {medication_2['rationale']}\n")

    print(f"3. Medication: {medication_3['name']} ({medication_3['class']})")
    print(f"   Rationale: {medication_3['rationale']}\n")

    print("-" * 80)

if __name__ == '__main__':
    recommend_htn_medications()