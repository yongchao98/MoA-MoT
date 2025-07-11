def recommend_hypertension_medications():
    """
    Analyzes the patient case and recommends a 3-drug regimen for resistant hypertension.
    """
    # Patient Profile
    patient_ethnicity = "African American"
    patient_a1c = 6.7
    patient_heart_rate = 91

    # Drug Recommendations and Rationale
    recommendations = [
        {
            "drug": "Amlodipine (Calcium Channel Blocker)",
            "rationale": (
                f"This is a first-line agent for treating hypertension in {patient_ethnicity} patients "
                "according to major guidelines. It is effective for blood pressure lowering and is not on the patient's exclusion list."
            )
        },
        {
            "drug": "Lisinopril (ACE Inhibitor)",
            "rationale": (
                f"Given the patient's A1C of {patient_a1c}%, indicating type 2 diabetes, an ACE inhibitor "
                "is recommended for its kidney-protective benefits in patients with both hypertension and diabetes."
            )
        },
        {
            "drug": "Metoprolol Succinate (Beta-Blocker)",
            "rationale": (
                f"The patient has an elevated heart rate of {patient_heart_rate} bpm (tachycardia). "
                "A cardioselective beta-blocker like metoprolol will help lower both her blood pressure and heart rate, "
                "addressing two issues with one medication."
            )
        }
    ]

    print("Based on the patient's profile, here are 3 recommended medications to maximize her hypertension treatment:\n")

    for i, rec in enumerate(recommendations, 1):
        print(f"Recommendation {i}: {rec['drug']}")
        print(f"Rationale: {rec['rationale']}\n")

if __name__ == "__main__":
    recommend_hypertension_medications()