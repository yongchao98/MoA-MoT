def recommend_htn_medications():
    """
    This function recommends a 3-drug regimen for the patient's resistant hypertension
    based on their clinical profile and medication exclusions.
    """
    # Rationale for choices:
    # 1. Losartan (ARB): First-line for patients with diabetes to protect the kidneys.
    # 2. Amlodipine (CCB): First-line for African American patients and synergistic with an ARB.
    # 3. Spironolactone (MRA): Recommended add-on therapy for resistant hypertension, especially when
    #    other diuretics are not an option. It is not on the patient's exclusion list.
    
    recommendation = {
        1: "Losartan",
        2: "Amlodipine",
        3: "Spironolactone"
    }

    print("Based on the patient's profile, here are 3 recommended medications to maximize hypertension treatment:")
    for number, drug in recommendation.items():
        print(f"{number}. {drug}")

if __name__ == '__main__':
    recommend_htn_medications()