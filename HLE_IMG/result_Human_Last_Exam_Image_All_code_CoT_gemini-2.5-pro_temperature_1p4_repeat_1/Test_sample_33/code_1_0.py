def diagnose_ecg(rhythm, rate, qrs_width, qrs_morphology):
    """
    This function simulates the diagnostic process for a wide-complex tachycardia.
    """
    print("ECG Analysis Steps:")
    print(f"1. Rhythm Analysis: The rhythm is {rhythm}.")
    print(f"2. Rate Analysis: The ventricular rate is {rate}.")
    print(f"3. QRS Analysis: The QRS width is {qrs_width} and its morphology is {qrs_morphology}.")

    # Diagnostic Logic
    diagnosis = "Undetermined"
    if rhythm == "irregularly irregular":
        print("\nReasoning: An 'irregularly irregular' rhythm strongly suggests Atrial Fibrillation.")
        if qrs_width == "wide":
            print("The wide QRS complexes indicate an abnormal ventricular conduction pathway.")
            if qrs_morphology == "variable and bizarre":
                print("The combination of AF, a very fast rate, and variable/bizarre QRS morphology is the classic presentation of conduction over an accessory pathway.")
                diagnosis = "Pre-excited Atrial Fibrillation"
            else:
                diagnosis = "Atrial Fibrillation with Aberrancy (e.g., Bundle Branch Block)"
        else:
            diagnosis = "Atrial Fibrillation with normal conduction"
    elif rhythm == "regular":
        diagnosis = "Could be SVT with aberrancy or Ventricular Tachycardia."
    
    print("\n--- Conclusion ---")
    print(f"The findings lead to a diagnosis of: {diagnosis}")
    return diagnosis

# Input features from the provided ECG
ecg_features = {
    "rhythm": "irregularly irregular",
    "rate": "very fast (>150 bpm)",
    "qrs_width": "wide (>0.12s)",
    "qrs_morphology": "variable and bizarre"
}

# Run the diagnostic function
final_diagnosis = diagnose_ecg(**ecg_features)