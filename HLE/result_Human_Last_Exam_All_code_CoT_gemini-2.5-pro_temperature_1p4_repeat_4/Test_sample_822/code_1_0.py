import sys

def diagnose_case():
    """
    Analyzes the clinical case to determine the most likely disease
    causing the patient's final illness.
    """

    # Define key data from the case
    patient_age = 62
    smoking_history_pack_years = 20

    print("Step-by-step Diagnosis based on the Clinical Case:")
    print("-" * 50)

    # Step 1: Analyze Patient Profile and Initial Presentation
    print("\n[Step 1] Analyzing Patient Profile and History:")
    print(f"- A {patient_age}-year-old male with a {smoking_history_pack_years}-pack-year smoking history.")
    print("- Occupation in shipbuilding suggests a risk of asbestos exposure.")
    print("- Initial symptoms (polyarthritis, fatigue) and pulmonary nodules led to steroid treatment.")
    print("- Conclusion: The patient was made immunocompromised by the steroid therapy.")

    # Step 2: Analyze the Acute Terminal Illness
    print("\n[Step 2] Analyzing the Acute Terminal Illness:")
    print("- The patient developed a severe infection with fever, productive cough, shortness of breath, and cutaneous lesions.")
    print("- This represents a systemic infection originating from the lungs and disseminating to the skin.")
    print("- The patient's confusion also suggests possible Central Nervous System (CNS) involvement.")
    print("- CRITICAL CLUE: The prescribed Aminoglycoside therapy was ineffective.")

    # Step 3: Synthesize and Conclude
    print("\n[Step 3] Synthesizing the Clues for a Final Diagnosis:")
    print("- The patient is an immunocompromised host.")
    print("- The infection affects multiple systems (pulmonary, cutaneous, likely CNS).")
    print("- The pathogen is resistant to aminoglycosides.")
    print("- This specific combination of factors is the classic presentation for an opportunistic infection caused by the bacterium *Nocardia*.")

    # Final Diagnosis
    print("\n[Final Conclusion]")
    print("The disease that best explains the patient's terminal illness—the severe, disseminated infection leading to septic shock—is Nocardiosis.")


if __name__ == '__main__':
    diagnose_case()