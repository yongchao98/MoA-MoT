def analyze_ecg():
    """
    Analyzes the provided ECG findings to determine the diagnosis.
    """
    print("ECG Analysis Steps:")
    
    print("\nStep 1: Assess Rhythm and Rate")
    print(" - Rhythm: The R-R intervals are variable, indicating an irregular rhythm.")
    print(" - Rate: The rate is rapid, approximately 150-200 beats per minute. This is a tachycardia.")
    print(" - QRS: The QRS complexes are wide (duration > 0.12 seconds).")
    print(" - Initial Impression: Irregular Wide Complex Tachycardia.")

    print("\nStep 2: Look for features of Ventricular Tachycardia (VT)")
    print(" - In any Wide Complex Tachycardia, VT should be assumed until proven otherwise.")
    print(" - Feature 1 (Lead aVR): There is a prominent, initial R wave in lead aVR. This is a highly specific sign for VT.")
    print(" - Feature 2 (Precordial Leads V1-V6): All QRS complexes from V1 through V6 are predominantly negative. This is called 'negative concordance' and is also highly specific for VT.")

    print("\nStep 3: Evaluate other possibilities")
    print(" - Atrial Fibrillation with Aberrancy: While it causes an irregular WCT, it would not typically show a dominant R in aVR or precordial concordance.")
    print(" - Pre-excited Atrial Fibrillation: Also an irregular WCT, but the specific VT criteria seen here make it less likely.")
    print(" - SVT with Aberrancy: This is typically a regular rhythm.")

    print("\nStep 4: Conclusion")
    print(" - The combination of a wide complex tachycardia with a dominant R wave in aVR and negative precordial concordance strongly indicates a ventricular origin.")
    print(" - Therefore, the diagnosis is Ventricular Tachycardia.")

analyze_ecg()