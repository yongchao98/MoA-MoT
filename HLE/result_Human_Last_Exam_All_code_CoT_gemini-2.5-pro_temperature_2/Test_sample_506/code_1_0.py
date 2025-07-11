def recommend_htn_medications():
    """
    Recommends a 3-drug regimen for the patient's resistant hypertension
    based on their profile and medication restrictions.
    """
    
    # Rationale:
    # 1. Lisinopril: ACE inhibitor, indicated for hypertension with comorbid diabetes for renal protection.
    # 2. Amlodipine: Dihydropyridine Calcium Channel Blocker, a first-line agent for African American patients.
    # 3. Carvedilol: A beta-blocker with alpha-blocking properties. Chosen as a third agent
    #    to provide an additional mechanism of action for blood pressure control and to help
    #    manage the patient's elevated heart rate.
    
    medication_recommendations = {
        "1. Lisinopril": "An ACE inhibitor to protect the kidneys due to diabetes and treat high blood pressure.",
        "2. Amlodipine": "A calcium channel blocker, which is a first-line and highly effective treatment for hypertension in African American patients.",
        "3. Carvedilol": "A beta-blocker with alpha-blocking properties to provide a third mechanism for blood pressure reduction and help control the patient's elevated heart rate."
    }

    print("Based on the patient's profile and medication restrictions, here is a recommended 3-drug regimen to maximize hypertension treatment:")
    for med, reason in medication_recommendations.items():
        print(f"\n{med}")
        print(f"   - Rationale: {reason}")

if __name__ == "__main__":
    recommend_htn_medications()