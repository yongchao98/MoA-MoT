def solve_medical_case():
    """
    Analyzes a clinical case to determine the most likely diagnosis.
    """

    # 1. Key findings extracted from the clinical case description
    findings = {
        "History": "62-year-old man, 20-pack-year smoker, ship building (asbestos exposure risk).",
        "Chronic Symptoms": "Polyarthritis (wrists, ankles, elbows), fatigue, confusion, bruising, dysphagia.",
        "Radiology": "Chest X-ray shows multiple pulmonary nodules.",
        "Treatment & Course": "Started on steroids, suggesting an inflammatory/autoimmune condition.",
        "Acute Illness": "After travel, developed fever, productive cough, cutaneous lesions, and fatal septic shock.",
        "Immunocompromise": "Steroid use created an immunocompromised state, leading to an overwhelming infection that was resistant to standard aminoglycoside therapy."
    }

    # 2. Differential Diagnoses with their explanatory power for the findings
    diagnoses = {
        "Granulomatosis with Polyangiitis (GPA)": 0,
        "Lung Cancer w/ Paraneoplastic Syndrome": 0,
        "Rheumatoid Arthritis w/ lung nodules": 0,
        "Sarcoidosis": 0
    }

    print("Analyzing the clinical case step-by-step:\n")

    # --- Scoring Logic ---

    # Finding 1: Pulmonary Nodules + Polyarthritis (joint pain/swelling)
    print("Step 1: Evaluating the combination of lung nodules and systemic joint inflammation.")
    # GPA is a classic cause of both lung nodules and inflammatory arthritis.
    diagnoses["Granulomatosis with Polyangiitis (GPA)"] += 2
    # Lung cancer can cause nodules, and paraneoplastic syndromes can cause joint pain.
    diagnoses["Lung Cancer w/ Paraneoplastic Syndrome"] += 1
    # RA can cause both but is a less complete explanation for the other symptoms.
    diagnoses["Rheumatoid Arthritis w/ lung nodules"] += 1
    # Sarcoidosis also fits this part of the picture.
    diagnoses["Sarcoidosis"] += 1
    print("Result: GPA is a strong match for this combination.\n")

    # Finding 2: Multi-systemic Symptoms (confusion, bruising, dysphagia)
    print("Step 2: Considering the broad range of other systemic symptoms.")
    # Vasculitis, the hallmark of GPA, affects blood vessels throughout the body, providing a unified explanation.
    diagnoses["Granulomatosis with Polyangiitis (GPA)"] += 2
    # Paraneoplastic syndromes can be very diverse and can explain these.
    diagnoses["Lung Cancer w/ Paraneoplastic Syndrome"] += 1
    # These are less typical for RA or Sarcoidosis.
    diagnoses["Rheumatoid Arthritis w/ lung nodules"] += 0
    diagnoses["Sarcoidosis"] += 0
    print("Result: The vasculitis of GPA provides the most direct explanation for these widespread symptoms.\n")

    # Finding 3: Severe Immunocompromise and Fatal Sepsis after Steroid Use
    print("Step 3: Analyzing the final illness.")
    # GPA is an autoimmune disease. Both the disease and its treatment (steroids) cause severe immunosuppression, making the patient extremely vulnerable to opportunistic infections, which are often resistant to standard antibiotics.
    diagnoses["Granulomatosis with Polyangiitis (GPA)"] += 2
    # Cancer patients are also immunocompromised, and steroids would worsen this.
    diagnoses["Lung Cancer w/ Paraneoplastic Syndrome"] += 1
    print("Result: The course of fatal infection after steroid treatment is highly characteristic of an underlying systemic autoimmune disease like GPA.\n")

    # Finding 4: Risk Factors (Smoking, Asbestos)
    print("Step 4: Considering patient history.")
    # These are very strong risk factors for lung cancer.
    diagnoses["Lung Cancer w/ Paraneoplastic Syndrome"] += 2
    # They are not specifically associated with GPA, although smoking can worsen many lung conditions.
    diagnoses["Granulomatosis with Polyangiitis (GPA)"] += 0
    print("Result: Patient history strongly suggests a primary lung malignancy.\n")
    
    # --- Conclusion ---
    print("--- Final Diagnosis ---")
    print("Scores:")
    for disease, score in diagnoses.items():
        print(f"- {disease}: {score}")

    # Determine the highest scoring diagnosis
    best_diagnosis = max(diagnoses, key=diagnoses.get)

    print("\nConclusion:")
    print("While the patient's history of smoking and asbestos exposure are major risk factors for Lung Cancer,")
    print("the specific clinical presentation of pulmonary nodules combined with polyarthritis, multi-system vasculitic symptoms (confusion, bruising), and extreme vulnerability to infection points more strongly toward Granulomatosis with Polyangiitis (GPA), also known as Wegener's Granulomatosis.")
    print(f"\nThe most likely diagnosis based on the overall clinical picture is: {best_diagnosis}")
    
    # Required final answer format
    return best_diagnosis

final_answer = solve_medical_case()
# <<<answer>>> is a placeholder for the final answer format requested in the prompt.
print(f"\n<<<Granulomatosis with Polyangiitis (Wegener's Granulomatosis)>>>")
