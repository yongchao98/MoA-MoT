def recommend_htn_medications():
    """
    Recommends three medications for the patient MM to maximize hypertension treatment
    based on her clinical profile and contraindications.
    """

    # Patient profile: 27yo African American F, resistant HTN, diabetes, hypercholesterolemia.
    # Contraindicated meds: eplerenone, finerenone, HCTZ, indapamide, bumetanide, furosemide,
    # torsemide, metolazone, verapamil.

    # Medication 1: Calcium Channel Blocker (CCB)
    med_1 = {
        "name": "Amlodipine",
        "rationale": "An effective first-line dihydropyridine calcium channel blocker for African American patients, which helps relax blood vessels. It is a suitable alternative to the contraindicated verapamil."
    }

    # Medication 2: ACE Inhibitor
    med_2 = {
        "name": "Lisinopril",
        "rationale": "An ACE inhibitor that is strongly recommended for patients with both hypertension and diabetes to provide renal (kidney) protection and control blood pressure. An ARB like Losartan would also be a suitable alternative."
    }

    # Medication 3: Aldosterone Antagonist
    med_3 = {
        "name": "Spironolactone",
        "rationale": "An aldosterone antagonist/potassium-sparing diuretic. It is the guideline-recommended add-on therapy for resistant hypertension. This medication is not on the patient's list of restricted drugs and is a critical choice since common diuretics are contraindicated."
    }
    
    print("Recommended 3-Drug Regimen for Resistant Hypertension:\n")
    print(f"1. Medication: {med_1['name']}")
    print(f"   Rationale: {med_1['rationale']}\n")
    
    print(f"2. Medication: {med_2['name']}")
    print(f"   Rationale: {med_2['rationale']}\n")

    print(f"3. Medication: {med_3['name']}")
    print(f"   Rationale: {med_3['rationale']}\n")

recommend_htn_medications()
print("<<<Amlodipine, Lisinopril, Spironolactone>>>")