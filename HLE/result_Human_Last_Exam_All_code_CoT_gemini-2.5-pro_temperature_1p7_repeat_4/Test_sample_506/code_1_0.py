def recommend_htn_medications():
    """
    Recommends a 3-drug regimen for the patient MM to treat resistant hypertension.
    """
    # Patient MM's relevant data points
    age = 27
    race = "African American"
    current_bp_systolic = "145-153"
    current_bp_diastolic = "85-89"
    goal_bp_systolic = 130
    goal_bp_diastolic = 90
    a1c = 6.7
    tsh = 4.5

    # Recommended Medications
    medications = {
        "Amlodipine": f"A dihydropyridine calcium channel blocker (CCB). This is a first-line agent for treating hypertension in {race} patients. It will help lower her BP from its current range of {current_bp_systolic}/{current_bp_diastolic} mmHg towards the goal of <{goal_bp_systolic}/{goal_bp_diastolic} mmHg.",
        "Losartan": f"An angiotensin II receptor blocker (ARB). This class is highly recommended for patients with hypertension and diabetes (patient's A1C is {a1c}%) as it provides crucial kidney protection.",
        "Spironolactone": "A mineralocorticoid receptor antagonist (MRA) that also functions as a potassium-sparing diuretic. It is a guideline-recommended medication for resistant hypertension, which is often difficult to control. This is a key choice given the patient is unable to take other common diuretic classes."
    }

    print(f"Based on the case of MM, a {age}-year-old female with resistant hypertension, the following 3 medications are recommended:")
    print("---")
    
    count = 1
    for med, rationale in medications.items():
        print(f"{count}. Medication: {med}")
        print(f"   Rationale: {rationale}\n")
        count += 1
        
if __name__ == '__main__':
    recommend_htn_medications()
    # The final answer containing the 3 recommended medications.
    print("<<<Amlodipine, Losartan, Spironolactone>>>")
