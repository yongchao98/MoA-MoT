def provide_diagnosis():
    """
    This function analyzes the provided clinical case and prints the most likely diagnosis and its justification.
    """
    # Patient Data from the Case
    patient_age = 62
    smoking_history = 20 # in pack-years

    diagnosis = "Granulomatosis with Polyangiitis (GPA), formerly known as Wegener's Granulomatosis"

    print(f"The most likely disease is: {diagnosis}")
    print("\n--- Diagnostic Justification ---")
    print(f"The patient, a {patient_age}-year-old man with a {smoking_history}-pack-year smoking history, presented a clinical picture that strongly points to GPA for several reasons:")
    print("\n1. Multisystem Presentation: The initial symptoms of joint pain (polyarthritis), fatigue, and later neurological changes (dizziness, confusion), bruising, and swallowing difficulty indicate a systemic disease affecting multiple organs, which is characteristic of a vasculitis like GPA.")
    print("\n2. Hallmark Lung Findings: The chest X-ray showed multiple pulmonary nodules. This is a classic and highly specific finding for GPA, representing areas of granulomatous inflammation in the lungs.")
    print("\n3. Complication of Immunosuppression: Patients with GPA are often treated with steroids. This treatment suppresses the immune system, making them vulnerable to severe opportunistic infections. The final illness after traveling to Africa (fever, productive cough, cutaneous lesions) fits the profile of an aggressive infection that took hold due to his compromised immune state.")
    print("\n4. Terminal Event: The patient's death from septic shock due to an infection that did not respond to standard aminoglycoside therapy is a known and unfortunate outcome in immunocompromised individuals. The pathogen was likely resistant or not a typical bacterium covered by that antibiotic class.")

provide_diagnosis()