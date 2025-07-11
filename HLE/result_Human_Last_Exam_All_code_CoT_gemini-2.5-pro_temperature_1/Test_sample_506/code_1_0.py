def recommend_hypertension_medications():
    """
    This function analyzes the patient's case and recommends a 3-drug regimen
    for resistant hypertension based on established clinical guidelines.
    """

    # Patient Profile
    patient_info = {
        "Age": 27,
        "Race": "African American",
        "Conditions": ["Resistant Hypertension Stage II", "Type 2 Diabetes (A1C 6.7%)", "Hypercholesterolemia", "High Heart Rate (91)"],
        "BP_Goal": "<130/90 mmHg"
    }

    # Medications the patient cannot take
    exclusion_list = [
        "eplerenone", "finerenone", "hydrochlorothiazide", "indapamide",
        "bumetanide", "furosemide", "torsemide", "metolazone", "verapamil"
    ]

    # Recommendation Logic
    recommendations = []
    reasoning = []

    # Medication 1: ARB for Diabetes and HTN
    med1 = "Losartan"
    reason1 = ("1. Losartan (Angiotensin II Receptor Blocker - ARB): This is a first-line agent for hypertension. "
               "It is strongly recommended for patients with Type 2 Diabetes as it provides crucial kidney protection.")
    recommendations.append(med1)
    reasoning.append(reason1)

    # Medication 2: CCB for Race and Synergy
    med2 = "Amlodipine"
    reason2 = ("2. Amlodipine (Dihydropyridine Calcium Channel Blocker - CCB): This is a first-line agent and is "
               "particularly effective for African American patients. It works synergistically with an ARB like Losartan to lower blood pressure.")
    recommendations.append(med2)
    reasoning.append(reason2)

    # Medication 3: Aldosterone Antagonist for Resistant HTN
    med3 = "Spironolactone"
    reason3 = ("3. Spironolactone (Aldosterone Antagonist): The patient has resistant hypertension, which typically requires a diuretic. "
               "Since common diuretics are on the exclusion list, Spironolactone is the ideal choice. It is a guideline-recommended "
               "agent for treating resistant hypertension and is not on the patient's exclusion list.")
    recommendations.append(med3)
    reasoning.append(reason3)

    # Print the final recommendations and reasoning
    print("Based on the patient's profile, the following 3 medications are recommended to maximize her hypertension treatment:\n")
    for reason in reasoning:
        print(reason + "\n")

    final_recommendation_str = ", ".join(recommendations)
    print("Final Recommended Regimen: " + final_recommendation_str)


if __name__ == "__main__":
    recommend_hypertension_medications()
<<<Losartan, Amlodipine, Spironolactone>>>