def recommend_hypertension_medications():
    """
    Analyzes a patient case and recommends a 3-drug regimen for resistant hypertension.
    """
    # Patient Data
    patient_info = {
        "age": 27,
        "ethnicity": "African American",
        "bp": "145-153/85-89",
        "hr": 91,
        "a1c_percent": 6.7,
        "is_diabetic": True, # Based on A1C > 6.5%
        "contraindications": [
            "eplerenone", "finerenone", "hydrochlorothiazide", "indapamide",
            "bumetanide", "furosemide", "torsemide", "metolazone", "verapamil"
        ]
    }

    recommendations = []

    # Recommendation 1: Calcium Channel Blocker (CCB)
    # Rationale: First-line for African American patients per JNC8/ACC-AHA guidelines.
    # A dihydropyridine CCB is appropriate as the non-dihydropyridine verapamil is excluded.
    med1_justification = (
        "1. Amlodipine (Calcium Channel Blocker): This is a recommended first-line agent for a patient "
        f"of {patient_info['ethnicity']} ethnicity with a blood pressure of {patient_info['bp']} mmHg. "
        "It is a potent vasodilator."
    )
    recommendations.append(med1_justification)

    # Recommendation 2: Angiotensin Receptor Blocker (ARB)
    # Rationale: Compelling indication for an ARB or ACEi due to diabetes for renal protection.
    med2_justification = (
        "2. Losartan (Angiotensin Receptor Blocker): This medication is indicated for kidney protection "
        f"in patients with diabetes, which is present here given an A1C of {patient_info['a1c_percent']}%. "
        "Combining an ARB with a CCB is a very effective strategy."
    )
    recommendations.append(med2_justification)

    # Recommendation 3: Beta-Blocker
    # Rationale: Diuretics and MRAs are excluded. A beta-blocker is a logical third-line agent
    # for resistant hypertension, especially given the patient's elevated heart rate.
    # Carvedilol has added alpha-blocking (vasodilating) properties.
    med3_justification = (
        "3. Carvedilol (Beta-Blocker with alpha-blocking activity): Given the exclusions of all diuretic classes and MRAs, "
        "a beta-blocker is a strong choice for a third agent. Carvedilol will help manage her elevated "
        f"heart rate of {patient_info['hr']} BPM and adds another mechanism of action for blood pressure control."
    )
    recommendations.append(med3_justification)

    print("Based on the patient's profile, the following three medications are recommended to maximize hypertension treatment:\n")
    for rec in recommendations:
        print(rec)
        print("-" * 20)

# Execute the recommendation function
if __name__ == "__main__":
    recommend_hypertension_medications()