def recommend_htn_medications():
    """
    This function analyzes the patient's case and recommends a three-drug regimen
    for hypertension based on clinical guidelines and patient-specific factors.
    """
    
    # Patient Profile:
    # - 27 y/o African American female
    # - Resistant Stage II Hypertension (BP: 145-153/85-89, Goal: <130/90)
    # - Type 2 Diabetes (A1C: 6.7%)
    # - Tachycardia (HR: 91)
    # - Contraindicated meds: eplerenone, finerenone, hydrochlorothiazide, indapamide,
    #   bumetanide, furosemide, torsemide, metolazone, verapamil.
    
    # Medication selection rationale:
    # 1. CCB (Amlodipine): First-line for African American patients. Dihydropyridine CCBs are appropriate as non-DHP (Verapamil) is excluded.
    # 2. ACE Inhibitor (Lisinopril): Compelling indication due to comorbid diabetes for renal protection.
    # 3. Beta-Blocker (Carvedilol): Addresses tachycardia and adds a third mechanism for resistant HTN.
    
    med1 = "Amlodipine"
    med2 = "Lisinopril"
    med3 = "Carvedilol"
    
    print("Based on the patient's clinical profile, here are three recommended medications to maximize her hypertension treatment:")
    print("-" * 40)
    
    # Printing each recommended medication
    print(f"1. {med1}: A dihydropyridine calcium channel blocker. This is a first-line agent that works by relaxing blood vessels.")
    print(f"2. {med2}: An ACE inhibitor. This is recommended due to the patient's diabetes, providing both blood pressure control and kidney protection.")
    print(f"3. {med3}: A beta-blocker with alpha-blocking activity. This addresses the patient's high heart rate and adds another mechanism to control resistant hypertension.")
    print("-" * 40)

    # Final answer format as requested
    final_answer = f"{med1}, {med2}, {med3}"
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    recommend_htn_medications()