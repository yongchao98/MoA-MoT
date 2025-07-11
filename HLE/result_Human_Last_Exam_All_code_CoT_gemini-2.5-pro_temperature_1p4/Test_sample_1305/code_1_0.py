import collections

def evaluate_cardiac_case():
    """
    This function analyzes a clinical vignette to determine the most likely diagnosis.
    It codifies the patient's findings and compares them against a knowledge base
    of cardiac conditions to find the best match.
    """

    # Step 1: Define the key findings from the patient's case
    patient_findings = {
        'murmur_type': 'systolic_ejection',
        'murmur_location': 'left_upper_sternal_border',
        'inspiration_effect': 'increase',
        'ecg_hypertrophy': 'RVH',
        'ecg_axis': 'LAD'
    }

    # Step 2: Create a knowledge base of findings for each diagnosis
    # Note: 'any' means the finding can be variable or not a key differentiator.
    # A list (e.g., ['RAD', 'LAD']) means either finding is possible depending on the subtype.
    diagnoses = {
        'A. Ebstein anomaly': {
            'murmur_type': 'holosystolic', # Mismatch
            'murmur_location': 'left_lower_sternal_border', # Mismatch
            'inspiration_effect': 'increase',
            'ecg_hypertrophy': 'RVH',
            'ecg_axis': 'RAD' # Mismatch
        },
        'B. Patent ductus arteriosus': {
            'murmur_type': 'continuous', # Mismatch
            'murmur_location': 'left_infraclavicular', # Mismatch
            'inspiration_effect': 'any',
            'ecg_hypertrophy': 'LVH', # Mismatch
            'ecg_axis': 'any'
        },
        'C. Mitral valve prolapse': {
            'murmur_type': 'late_systolic', # Mismatch
            'murmur_location': 'apex', # Mismatch
            'inspiration_effect': 'any',
            'ecg_hypertrophy': 'none', # Mismatch
            'ecg_axis': 'any'
        },
        'D. Atrial septal defect': {
            'murmur_type': 'systolic_ejection',
            'murmur_location': 'left_upper_sternal_border',
            'inspiration_effect': 'increase',
            'ecg_hypertrophy': 'RVH',
            'ecg_axis': ['RAD', 'LAD'] # Match (LAD is classic for ostium primum ASD)
        },
        'E. Hypertrophic cardiomyopathy': {
            'murmur_type': 'systolic_ejection',
            'murmur_location': 'left_lower_sternal_border', # Mismatch
            'inspiration_effect': 'decrease', # Mismatch
            'ecg_hypertrophy': 'LVH', # Mismatch
            'ecg_axis': 'LAD'
        },
        'F. Tricuspid stenosis': {
            'murmur_type': 'diastolic', # Mismatch
            'murmur_location': 'left_lower_sternal_border', # Mismatch
            'inspiration_effect': 'increase',
            'ecg_hypertrophy': 'none', # Mismatch (RAE is more common)
            'ecg_axis': 'RAD' # Mismatch
        },
        'G. Ventricular septal defect': {
            'murmur_type': 'holosystolic', # Mismatch
            'murmur_location': 'left_lower_sternal_border', # Mismatch
            'inspiration_effect': 'any',
            'ecg_hypertrophy': 'LVH', # Mismatch
            'ecg_axis': 'any'
        }
    }

    # Step 3: Score each diagnosis against the patient's findings
    results = {}
    print("--- Diagnostic Scorecard ---")
    print(f"Patient Findings: {patient_findings}\n")

    for name, profile in diagnoses.items():
        score = 0
        matches = []
        for key, value in patient_findings.items():
            profile_value = profile.get(key)
            # Check for a match. Handles list of possible values (like for ASD axis).
            if isinstance(profile_value, list):
                if value in profile_value:
                    score += 1
                    matches.append(key)
            elif profile_value == value:
                score += 1
                matches.append(key)
        results[name] = score
        print(f"Diagnosis: {name}")
        print(f"  - Matching features: {matches}")
        print(f"  - Match Score: {score} out of {len(patient_findings)}\n")

    # Step 4: Identify and print the best match
    best_match = max(results, key=results.get)
    print("--- Conclusion ---")
    print("The analysis shows that Atrial Septal Defect (D) is the only diagnosis that can account for all the patient's key findings.")
    print("Specifically, an ostium primum type of ASD perfectly explains the combination of:")
    print("1. A systolic ejection murmur at the left upper sternal border (from relative pulmonic stenosis).")
    print("2. Signs of right heart volume overload (RVH on ECG).")
    print("3. The highly specific ECG finding of Left Axis Deviation (LAD) combined with RVH.")
    print(f"\nThe most likely diagnosis is: {best_match}")

evaluate_cardiac_case()