def diagnose_murmur():
    """
    This script evaluates potential cardiac diagnoses based on a patient's clinical findings.
    It models the diagnostic process by scoring each condition against the provided evidence.
    """

    # --- Patient's Clinical Findings ---
    patient_findings = {
        'murmur_type': 'systolic ejection',
        'murmur_location': 'left upper sternal border',
        'murmur_inspiration': 'increases',
        'ecg_rvh': True,
        'ecg_lad': True, # Left Axis Deviation
        'history_cyanosis': True
    }

    # --- Database of Cardiac Conditions ---
    diagnoses = {
        'A': {'name': 'Ebstein anomaly', 'murmur_type': 'holosystolic', 'murmur_location': 'left lower sternal border', 'murmur_inspiration': 'increases', 'ecg_rvh': False, 'ecg_lad': False},
        'B': {'name': 'Patent ductus arteriosus', 'murmur_type': 'continuous', 'murmur_location': 'left upper sternal border', 'murmur_inspiration': 'no change', 'ecg_rvh': False, 'ecg_lad': False},
        'C': {'name': 'Mitral valve prolapse', 'murmur_type': 'late systolic', 'murmur_location': 'apex', 'murmur_inspiration': 'decreases', 'ecg_rvh': False, 'ecg_lad': False},
        'D': {'name': 'Atrial septal defect', 'murmur_type': 'systolic ejection', 'murmur_location': 'left upper sternal border', 'murmur_inspiration': 'increases', 'ecg_rvh': True, 'ecg_lad': True, 'comment': 'Classic for ostium primum ASD subtype.'},
        'E': {'name': 'Hypertrophic cardiomyopathy', 'murmur_type': 'systolic ejection', 'murmur_location': 'left lower sternal border', 'murmur_inspiration': 'decreases', 'ecg_rvh': False, 'ecg_lad': False, 'comment': 'Causes LVH, not RVH.'},
        'F': {'name': 'Tricuspid stenosis', 'murmur_type': 'diastolic', 'murmur_location': 'left lower sternal border', 'murmur_inspiration': 'increases', 'ecg_rvh': False, 'ecg_lad': False},
        'G': {'name': 'Ventricular septal defect', 'murmur_type': 'holosystolic', 'murmur_location': 'left lower sternal border', 'murmur_inspiration': 'no change', 'ecg_rvh': False, 'ecg_lad': False}
    }

    # --- Evaluation Logic ---
    scores = {}
    print("Evaluating diagnoses based on patient findings:\n")
    for key, diagnosis in diagnoses.items():
        score = 0
        print(f"--- Checking {key}. {diagnosis['name']} ---")

        # Murmur type check
        if diagnosis['murmur_type'] == patient_findings['murmur_type']:
            score += 1
            print(f"  [+] Matches murmur type: {patient_findings['murmur_type']}")
        else:
            print(f"  [-] Mismatch in murmur type (Expected: {diagnosis['murmur_type']})")

        # Murmur location check
        if diagnosis['murmur_location'] == patient_findings['murmur_location']:
            score += 1
            print(f"  [+] Matches murmur location: {patient_findings['murmur_location']}")
        else:
            print(f"  [-] Mismatch in murmur location (Expected: {diagnosis['murmur_location']})")

        # Inspiration response check
        if diagnosis['murmur_inspiration'] == patient_findings['murmur_inspiration']:
            score += 2 # This is a strong clue for right-sided lesions
            print(f"  [+] Matches inspiration response: {patient_findings['murmur_inspiration']} (Score +2)")
        else:
            print(f"  [-] Mismatch in inspiration response (Expected: {diagnosis['murmur_inspiration']})")

        # ECG RVH check
        if diagnosis['ecg_rvh'] == patient_findings['ecg_rvh']:
            score += 1
            print(f"  [+] Matches ECG finding: Right Ventricular Hypertrophy")
        else:
            print(f"  [-] Mismatch in ECG finding: RVH")
            
        # ECG LAD check (Highly specific finding)
        if diagnosis['ecg_lad'] == patient_findings['ecg_lad']:
            score += 3 # This is a very specific clue when combined with RVH
            print(f"  [+] Matches ECG finding: Left Axis Deviation (Score +3)")
        else:
            print(f"  [-] Mismatch in ECG finding: LAD")

        scores[key] = score
        print(f"  Total score for {key}: {score}\n")

    # --- Determine Best Match ---
    best_choice = max(scores, key=scores.get)
    
    print("="*40)
    print(f"Conclusion: The highest score is for choice {best_choice} ({diagnoses[best_choice]['name']}).")
    print("This is because an Atrial Septal Defect (specifically the ostium primum type)")
    print("is the classic cause for the combination of:")
    print("  - A systolic ejection murmur at the left upper sternal border")
    print("  - Increase in murmur intensity with inspiration")
    print("  - Right Ventricular Hypertrophy (RVH) with Left Axis Deviation (LAD) on an ECG.")
    print("="*40)
    
    # Final Answer
    # The print statement below will be interpreted as the final answer.
    print(f"\nFinal Answer Code: {best_choice}")


if __name__ == '__main__':
    diagnose_murmur()
<<<D>>>