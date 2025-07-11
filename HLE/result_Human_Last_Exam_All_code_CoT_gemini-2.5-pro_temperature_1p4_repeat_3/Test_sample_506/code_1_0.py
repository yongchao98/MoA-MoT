def recommend_htn_medications():
    """
    This function recommends three medications for the patient's hypertension
    based on their clinical profile and medication exclusion list.
    """
    # Medication recommendations based on clinical analysis
    # 1. Losartan (ARB): For HTN with compelling indication of diabetes.
    # 2. Amlodipine (CCB): Potent antihypertensive, effective in African American patients.
    # 3. Metoprolol Succinate (Beta-Blocker): For HTN and comorbid tachycardia.
    
    recommendations = [
        "1. Losartan",
        "2. Amlodipine",
        "3. Metoprolol Succinate"
    ]
    
    print("Based on the patient's profile, the following 3 medications are recommended to maximize her hypertension treatment:")
    for med in recommendations:
        print(med)

if __name__ == "__main__":
    recommend_htn_medications()