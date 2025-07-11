def diagnose_patient_case():
    """
    Analyzes a clinical case study to determine the most likely diagnosis.
    """
    # 1. Patient Data
    patient_profile = {
        "age": 62,
        "history": "20-pack-year smoking",
        "occupation": "ship building"
    }

    symptoms_and_findings = {
        "initial": "fatigue, polyarthritis (wrists, ankles, elbows)",
        "progressive": "dizziness, confusion, bruising, dysphagia, loss of appetite, shortness of breath",
        "imaging": "multiple pulmonary nodules on Chest X-ray",
        "acute_event": "fever, productive cough (green sputum), cutaneous lesions after travel",
        "outcome": "died from septic shock despite aminoglycoside therapy"
    }

    # 2. Diagnostic Reasoning
    print("Step 1: Analyzing patient risk factors.")
    print(f"The patient is a {patient_profile['age']}-year-old male with a {patient_profile['history']} history and occupational exposure from {patient_profile['occupation']}.")
    print("This profile raises suspicion for lung cancer, asbestosis, and autoimmune conditions.\n")

    print("Step 2: Evaluating the constellation of symptoms.")
    print("The key features are the combination of systemic inflammation (polyarthritis), lower respiratory tract disease (multiple pulmonary nodules), and multi-organ involvement (neurologic, cutaneous, constitutional).\n")

    print("Step 3: Considering the differential diagnosis.")
    print("- Lung Cancer: Possible, but less likely to cause this severe, inflammatory polyarthritis.")
    print("- Rheumatoid Arthritis: Explains arthritis and can cause lung nodules, but the bruising, confusion, and rapid decline suggest a more aggressive vasculitic process.")
    print("- Infection (e.g., TB, Fungal): Unlikely as a primary cause given the chronic initial phase, though a secondary infection was the terminal event.\n")

    print("Step 4: Synthesizing the findings to reach a conclusion.")
    print("The clinical picture strongly points to a systemic small-vessel vasculitis. Granulomatosis with Polyangiitis (GPA) classically presents with a triad of upper respiratory, lower respiratory, and kidney involvement, but its presentation is variable.")
    print("In this case:")
    print(" - Lower respiratory involvement is clear (pulmonary nodules, shortness of breath).")
    print(" - Systemic inflammation is clear (polyarthritis, fatigue).")
    print(" - Multi-organ vasculitis explains the bruising (skin), confusion (CNS), and possibly dysphagia.")
    print(" - Immunosuppression from the disease and steroid treatment made the patient vulnerable to the fatal infection that caused septic shock.\n")

    # 5. Final Diagnosis
    final_diagnosis = "Granulomatosis with Polyangiitis (GPA)"
    print("---CONCLUSION---")
    print(f"The patient's presentation is most consistent with a diagnosis of: {final_diagnosis}")

diagnose_patient_case()