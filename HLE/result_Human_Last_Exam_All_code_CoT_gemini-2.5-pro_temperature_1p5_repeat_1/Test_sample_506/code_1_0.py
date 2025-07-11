def recommend_htn_medications():
    """
    Recommends a 3-drug regimen for a patient with resistant hypertension,
    based on their clinical profile and contraindications.
    """

    recommendations = [
        {
            "medication_name": "Amlodipine",
            "drug_class": "Dihydropyridine Calcium Channel Blocker (CCB)",
            "rationale": "First-line agent for African American patients per guidelines. Since most diuretics are contraindicated for this patient, a CCB is an excellent foundational therapy for her hypertension."
        },
        {
            "medication_name": "Losartan",
            "drug_class": "Angiotensin II Receptor Blocker (ARB)",
            "rationale": "Strongly recommended for patients with both hypertension and diabetes (A1C > 6.5%) to provide kidney protection. It works synergistically with Amlodipine to lower blood pressure."
        },
        {
            "medication_name": "Carvedilol",
            "drug_class": "Mixed Alpha/Beta-Adrenergic Blocker",
            "rationale": "An effective agent for resistant hypertension, especially when diuretics cannot be used. It will also help manage the patient's elevated heart rate (tachycardia, HR=91)."
        }
    ]

    print("Recommended 3-Medication Regimen for HTN Maximization:\n")
    for i, med in enumerate(recommendations, 1):
        print(f"Medication {i}: {med['medication_name']}")
        print(f"  - Drug Class: {med['drug_class']}")
        print(f"  - Rationale: {med['rationale']}\n")

    # Construct the final answer string for the required format
    final_answer = ", ".join([med['medication_name'] for med in recommendations])
    print(f"<<<{final_answer}>>>")

# Execute the function
recommend_htn_medications()