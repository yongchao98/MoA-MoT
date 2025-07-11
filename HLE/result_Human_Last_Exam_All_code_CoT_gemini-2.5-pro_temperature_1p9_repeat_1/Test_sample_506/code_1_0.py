import sys

def recommend_htn_medications():
    """
    Analyzes the patient case and recommends three medications to maximize hypertension (HTN) treatment.
    """
    
    patient_profile = {
        "Ethnicity": "African American",
        "Condition": "Resistant Stage II HTN",
        "Comorbidities": "Diabetes (A1C=6.7%), Hypercholesterolemia",
        "Goal BP": "<130/90 mmHg",
        "Constraints": ["eplerenone", "finerenone", "hydrochlorothiazide", "indapamide", 
                        "bumetanide", "furosemide", "torsemide", "metolazone", "verapamil"]
    }

    # Medication Recommendation Logic
    medication_1 = {
        "name": "Losartan",
        "class": "Angiotensin II Receptor Blocker (ARB)",
        "rationale": "This is a first-line medication for hypertension. It is especially beneficial for MM as she has diabetes, and ARBs help protect the kidneys. It works by relaxing blood vessels."
    }
    
    medication_2 = {
        "name": "Amlodipine",
        "class": "Dihydropyridine Calcium Channel Blocker (CCB)",
        "rationale": "This is another first-line medication that is highly effective in African American patients. It also works by relaxing blood vessels to lower blood pressure."
    }

    medication_3 = {
        "name": "Spironolactone",
        "class": "Mineralocorticoid Receptor Antagonist (MRA) / Potassium-Sparing Diuretic",
        "rationale": "MM has resistant hypertension, and most common diuretics are not an option for her. Spironolactone is the recommended medication for resistant hypertension when a standard diuretic isn't sufficient or can't be used. It helps the body remove excess salt and water and is not on her list of restricted medications."
    }
    
    recommendations = [medication_1, medication_2, medication_3]

    print("Based on the patient's profile, here are 3 recommended medications to maximize her hypertension treatment:")
    print("-" * 40)
    
    for i, med in enumerate(recommendations, 1):
        print(f"Recommendation {i}: {med['name']} ({med['class']})")
        print(f"  - Rationale: {med['rationale']}\n")

    # Combine medication names for the final answer format
    final_answer_str = ", ".join([med['name'] for med in recommendations])
    # The line below is for the final answer extraction.
    # The file is configured to extract the answer from the '<<<' and '>>>' symbols.
    sys.stdout.write(f'<<<{final_answer_str}>>>')

if __name__ == '__main__':
    recommend_htn_medications()