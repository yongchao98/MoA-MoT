import sys

def solve_medical_mystery():
    """
    Analyzes a patient's clinical case to determine the most probable diagnosis.
    The function breaks down the patient's history, symptoms, and treatment course
    to synthesize a final conclusion.
    """
    
    # --- Step 1: Analyze Patient History and Risk Factors ---
    # The patient is older and has a history that predisposes to lung disease and immunosuppression.
    patient_profile = {
        "Age": 62,
        "Smoking": "20-pack-year history",
        "Occupation": "ship building (asbestos exposure)",
        "Immunosuppression": "Started on steroids"
    }
    
    print("Patient Profile and Key Risk Factors:")
    print(f"- Age: {patient_profile['Age']}")
    print(f"- Smoking History: {patient_profile['Smoking']}")
    print(f"- Occupational Exposure: {patient_profile['Occupation']}")
    print(f"- Iatrogenic Risk: {patient_profile['Immunosuppression']}\n")
    
    # --- Step 2: Evaluate the Clinical Presentation ---
    # The disease has features involving multiple organ systems, which is a hallmark of certain infections.
    clinical_features = {
        "Pulmonary": "Multiple pulmonary nodules, shortness of breath, productive cough",
        "Systemic": "Fatigue, fever, loss of appetite",
        "Central Nervous System (CNS)": "Dizziness, confusion",
        "Cutaneous": "Cutaneous lesions, bruising"
    }
    
    print("Clinical Presentation (Multi-systemic Involvement):")
    print(f"- Pulmonary: {clinical_features['Pulmonary']}")
    print(f"- Central Nervous System: {clinical_features['CNS']}")
    print(f"- Cutaneous (Skin): {clinical_features['Cutaneous']}\n")
    
    # --- Step 3: Analyze the Treatment Course and Outcome ---
    # The failure of a specific antibiotic therapy is a crucial diagnostic clue.
    treatment_failure = "Aminoglycoside therapy proved ineffective."
    
    print("Crucial Diagnostic Clue from Treatment:")
    print(f"- {treatment_failure}\n")

    # --- Step 4: Synthesize and Conclude ---
    # The diagnosis must explain all the elements:
    # 1. An opportunistic infection in an immunocompromised host (due to steroid use).
    # 2. An infection that classically presents with a triad of lung, skin, and brain involvement.
    # 3. An organism that is not susceptible to aminoglycoside antibiotics.
    # Nocardia species fit this description perfectly. Nocardiosis is a bacterial infection that
    # mimics fungal or TB infections. It typically affects immunocompromised individuals,
    # disseminates from the lungs to the brain and skin, and is treated with sulfonamides,
    # not typically aminoglycosides.
    
    final_conclusion = "The combination of an immunocompromised host, pulmonary nodules, and disseminated disease to the brain (confusion) and skin (lesions) that is refractory to aminoglycoside therapy is a classic presentation of Nocardiosis."
    
    print("Final Diagnostic Reasoning:")
    print(final_conclusion)
    
    # --- Final Answer ---
    final_disease = "Nocardiosis"
    print(f"\nTherefore, the most likely disease is: <<<{final_disease}>>>")

if __name__ == "__main__":
    solve_medical_mystery()