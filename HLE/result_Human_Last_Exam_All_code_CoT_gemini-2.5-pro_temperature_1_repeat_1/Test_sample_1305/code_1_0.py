def diagnose_murmur():
    """
    Analyzes a clinical vignette to determine the most likely cause of a heart murmur.
    """

    # Key findings from the clinical case
    patient_age = 31
    history = "Faint cyanosis as a child"
    murmur_type = "Systolic ejection murmur"
    murmur_location = "Left upper sternal border"
    respiration_effect = "Increases in intensity with inspiration"
    ecg_findings = ["Left axis deviation", "Right ventricular hypertrophy (RVH)"]

    print("Analyzing the clinical case step-by-step:")
    print("-----------------------------------------")
    print(f"1. Patient Presentation: A {patient_age}-year-old woman with shortness of breath, fatigue, and a history of childhood cyanosis.")
    print("\n2. Murmur Analysis:")
    print(f"   - The murmur is a '{murmur_type}' at the '{murmur_location}'.")
    print(f"   - It '{respiration_effect}'. This (Carvallo's sign) strongly indicates a right-sided heart pathology.")
    print("   - This murmur is consistent with either pulmonic stenosis or high flow across the pulmonic valve (e.g., from an Atrial Septal Defect).")

    print("\n3. ECG Analysis:")
    print(f"   - The ECG shows a combination of '{ecg_findings[0]}' and '{ecg_findings[1]}'.")
    print("   - This is a critical and paradoxical finding. RVH is expected from a right-sided volume/pressure overload, but this should cause RIGHT axis deviation.")
    print("   - The specific combination of Left Axis Deviation with Right Ventricular Hypertrophy is a classic finding for an ostium primum Atrial Septal Defect (ASD).")

    print("\n4. Synthesis and Conclusion:")
    print("   - An ostium primum ASD explains all findings: the right-sided flow murmur, the adult onset of symptoms, the history of cyanosis (from transient shunting), and most importantly, the unique ECG pattern.")
    print("   - Therefore, among the choices provided, Atrial Septal Defect is the most likely diagnosis.")
    print("-----------------------------------------")

# Run the diagnosis
diagnose_murmur()