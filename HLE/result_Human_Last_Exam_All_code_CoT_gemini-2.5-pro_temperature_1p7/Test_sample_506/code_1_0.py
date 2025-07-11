def recommend_htn_medications():
    """
    Recommends a 3-drug regimen for the patient MM's resistant hypertension
    based on her clinical profile and medication restrictions.
    """
    patient_profile = {
        "Ethnicity": "African American",
        "Key Conditions": ["Resistant Hypertension", "Type 2 Diabetes", "Hypercholesterolemia"],
        "Heart Rate": 91,
        "Contraindicated Drug Classes": ["Most Diuretics", "Specific MRAs", "Non-DHP CCBs"]
    }

    medications = {
        "Amlodipine": "A Dihydropyridine Calcium Channel Blocker (CCB). This is a first-line, potent antihypertensive agent that is particularly effective for African American patients.",
        "Lisinopril": "An ACE Inhibitor (ACEi). This is a crucial medication for a patient with both hypertension and diabetes, as it provides kidney protection and works on the renin-angiotensin system.",
        "Carvedilol": "A mixed Alpha/Beta-Blocker. This is an excellent choice for a third agent when diuretics are not an option. It will help control her elevated heart rate (91 bpm) and provides an additional mechanism for blood pressure reduction with a favorable metabolic profile."
    }

    print("Based on the patient's profile, here are three recommended medications to maximize her hypertension treatment:\n")
    
    med_number = 1
    for med, rationale in medications.items():
        print(f"Recommendation {med_number}: {med}")
        print(f"   - Rationale: {rationale}\n")
        med_number += 1

recommend_htn_medications()
<<<Amlodipine, Lisinopril, Carvedilol>>>