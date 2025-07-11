import sys

def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring potential conditions against the patient's findings.
    """
    # Step 1: Define the patient's key clinical findings.
    # The combination of Right Ventricular Hypertrophy (RVH) and Left Axis Deviation (LAD) is a crucial and uncommon finding.
    patient_findings = {
        "murmur_type": "systolic ejection",
        "murmur_location": "left upper sternal border",
        "murmur_dynamic": "increases with inspiration",
        "ecg_hypertrophy": "right ventricular hypertrophy",
        "ecg_axis": "left axis deviation",
        "history": "childhood cyanosis"
    }

    # Step 2: Create a knowledge base of findings for each potential diagnosis.
    diagnoses = {
        "A. Ebstein anomaly": {
            "murmur_type": "holosystolic (tricuspid regurgitation)", "murmur_location": "left lower sternal border",
            "murmur_dynamic": "increases with inspiration", "ecg_hypertrophy": "right atrial enlargement",
            "ecg_axis": "right axis deviation", "history": "childhood cyanosis"
        },
        "B. Patent ductus arteriosus": {
            "murmur_type": "continuous machine-like", "murmur_location": "left infraclavicular area",
            "murmur_dynamic": "no change with inspiration", "ecg_hypertrophy": "left ventricular hypertrophy",
            "ecg_axis": "normal or left axis deviation", "history": "none"
        },
        "C. Mitral valve prolapse": {
            "murmur_type": "mid-systolic click, late systolic murmur", "murmur_location": "apex",
            "murmur_dynamic": "decreases with inspiration", "ecg_hypertrophy": "none",
            "ecg_axis": "normal", "history": "none"
        },
        "D. Atrial septal defect": {
            # Note: A primum ASD specifically causes LAD + RVH. Other types cause RAD.
            # The murmur is a flow murmur across the pulmonic valve.
            "murmur_type": "systolic ejection", "murmur_location": "left upper sternal border",
            "murmur_dynamic": "increases with inspiration", "ecg_hypertrophy": "right ventricular hypertrophy",
            "ecg_axis": "left axis deviation", # This is classic for a primum ASD.
            "history": "childhood cyanosis" # Can occur with transient shunt reversal.
        },
        "E. Hypertrophic cardiomyopathy with obstruction": {
            "murmur_type": "systolic ejection", "murmur_location": "left lower sternal border",
            "murmur_dynamic": "decreases with inspiration", # Key differentiator
            "ecg_hypertrophy": "left ventricular hypertrophy", "ecg_axis": "variable, often left", "history": "none"
        },
        "F. Tricuspid stenosis": {
            "murmur_type": "mid-diastolic rumble", "murmur_location": "left lower sternal border",
            "murmur_dynamic": "increases with inspiration", "ecg_hypertrophy": "right atrial enlargement",
            "ecg_axis": "normal or right axis deviation", "history": "none"
        },
        "G. Ventricular septal defect": {
            "murmur_type": "holosystolic", "murmur_location": "left lower sternal border",
            "murmur_dynamic": "no change with inspiration", "ecg_hypertrophy": "left ventricular hypertrophy",
            "ecg_axis": "normal or left axis deviation", "history": "none"
        },
    }

    print("Analyzing patient findings against possible diagnoses:\n")
    
    best_match = ""
    highest_score = -1

    # Step 3 & 4: Compare patient findings to the knowledge base and print the analysis.
    for diagnosis, findings in diagnoses.items():
        score = 0
        print(f"--- Checking {diagnosis} ---")
        for key, value in patient_findings.items():
            if key in findings and findings[key] == value:
                score += 1
                print(f"  [+] Match: {key.replace('_', ' ').title()} = '{value}'")
            else:
                print(f"  [-] Mismatch: Patient has '{value}' for {key.replace('_', ' ').title()}, but expect '{findings.get(key, 'N/A')}'")
        
        print(f"  Score for {diagnosis}: {score} out of {len(patient_findings)}\n")

        if score > highest_score:
            highest_score = score
            best_match = diagnosis

    print("--- Conclusion ---")
    print(f"The highest score is {highest_score}.")
    print(f"The most likely cause of this woman's murmur is {best_match}.")
    print("This is because Atrial Septal Defect (specifically the primum type, part of an endocardial cushion defect) is the classic cause for the combined ECG findings of left axis deviation and right ventricular hypertrophy, along with a systolic ejection murmur at the left upper sternal border due to increased flow across the pulmonic valve.")
    # Extract the letter from the best match string, e.g., 'D. ...' -> 'D'
    final_answer_letter = best_match.split('.')[0]
    sys.stdout.flush() # Ensure all previous prints are displayed before the final answer
    print(f'<<<{final_answer_letter}>>>')


solve_medical_case()