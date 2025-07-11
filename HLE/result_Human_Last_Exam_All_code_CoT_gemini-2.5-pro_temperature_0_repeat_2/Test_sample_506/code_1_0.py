def recommend_htn_medications():
    """
    Recommends a 3-drug regimen for the patient's resistant hypertension
    based on their clinical profile and medication restrictions.
    """
    
    # Patient Profile: 27 y/o African American female with resistant HTN,
    # diabetes, and hypercholesterolemia.
    # Goal: Maximize HTN treatment with 3 medications.
    # Constraints: Cannot take eplerenone, finerenone, HCTZ, indapamide,
    # loop diuretics, metolazone, or verapamil.

    # Medication 1: ACE Inhibitor for RAAS blockade and renal protection in diabetes.
    medication_1 = "Lisinopril"
    reason_1 = "An ACE inhibitor to block the renin-angiotensin system. It is a first-line agent that is also recommended for patients with diabetes to protect kidney function."

    # Medication 2: Dihydropyridine Calcium Channel Blocker for potent vasodilation.
    medication_2 = "Amlodipine"
    reason_2 = "A dihydropyridine calcium channel blocker. This class of medication is particularly effective for treating hypertension in African American patients and works by relaxing blood vessels."

    # Medication 3: Mineralocorticoid Receptor Antagonist for resistant hypertension.
    medication_3 = "Spironolactone"
    reason_3 = "A mineralocorticoid receptor antagonist. This is a guideline-recommended add-on therapy for resistant hypertension. It works as a potassium-sparing diuretic to help reduce blood pressure when other medications are not sufficient."

    print("Recommended 3-Medication Regimen for Hypertension:")
    print("-" * 50)
    print(f"1. {medication_1}: {reason_1}")
    print(f"2. {medication_2}: {reason_2}")
    print(f"3. {medication_3}: {reason_3}")
    print("-" * 50)
    print("\nNote: This regimen targets three different pathways to effectively lower blood pressure. The patient's high cholesterol and undertreated hypothyroidism should also be addressed separately to improve overall cardiovascular health.")

recommend_htn_medications()