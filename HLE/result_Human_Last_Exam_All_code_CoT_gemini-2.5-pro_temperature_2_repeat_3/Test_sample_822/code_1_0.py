def diagnose_patient_case():
    """
    Analyzes a clinical case and presents a reasoned diagnosis.
    """
    case_summary = {
        "Age": 62,
        "Gender": "Male",
        "Risk Factors": ["20-pack-year smoking history", "Shipbuilding occupation (asbestos risk)"],
        "Initial Symptoms": ["Polyarthritis (wrists, ankles, elbows)", "Fatigue"],
        "Progressive Symptoms": [
            "Neurological (dizziness, confusion)",
            "Hematological (bruising)",
            "Swallowing difficulty (dysphagia)",
            "Pulmonary (shortness of breath)",
            "Constitutional (loss of appetite)"
        ],
        "Imaging": {
            "Chest X-ray": "Multiple pulmonary nodules"
        },
        "Clinical Course": [
            "Started on steroids and NSAIDs",
            "Developed acute infection after travel (fever, productive cough, cutaneous lesions)",
            "Antibiotic therapy (Aminoglycoside) was ineffective",
            "Died from septic shock"
        ]
    }

    print("Analyzing patient case based on provided data...\n")

    print("### Step 1: Evaluating Key Clinical Features")
    print(f"- Patient presents with multi-system disease at age {case_summary['Age']}.")
    print(f"- Key finding on X-ray: {case_summary['Imaging']['Chest X-ray']}.")
    print(f"- Initial symptoms include: {', '.join(case_summary['Initial Symptoms'])}.")
    print(f"- Disease progression involved multiple systems: {', '.join(case_summary['Progressive Symptoms'])}.\n")

    print("### Step 2: Correlating Symptoms with Potential Diagnoses")
    print("The combination of symptoms points towards a systemic inflammatory condition.")
    print("- Arthritis + Multiple Pulmonary Nodules: This combination is a hallmark of certain autoimmune diseases, especially systemic vasculitides.")
    print("- Other features like confusion, bruising, and difficulty swallowing suggest potential kidney, blood vessel, and upper airway inflammation.")
    print("- Smoking history and nodules could suggest lung cancer, but the prominent arthritis and other systemic signs make a rheumatological cause more likely.\n")

    print("### Step 3: Assessing the Final Clinical Outcome")
    print("The patient's death from septic shock after an acute infection is a critical clue.")
    print("- The underlying disease and the use of steroids would have made the patient immunocompromised.")
    print("- This explains the susceptibility to a severe, opportunistic infection that was refractory to standard antibiotics.\n")

    print("### Step 4: Conclusion")
    print("The most unifying diagnosis is Granulomatosis with Polyangiitis (GPA).")
    print("Reasoning:")
    print("1. GPA is a vasculitis that classically affects the upper airways, lower airways (lungs), and kidneys.")
    print("2. It explains the multiple pulmonary nodules, arthritis, and constitutional symptoms.")
    print("3. The progressive symptoms (confusion, bruising, dysphagia) align with multi-organ vasculitic inflammation.")
    print("4. Treatment for GPA involves immunosuppressants like steroids, which directly explains why the patient was highly vulnerable to a fatal opportunistic infection.")

# Execute the diagnostic analysis
diagnose_patient_case()