def recommend_htn_medications():
    """
    Analyzes the patient case and recommends a 3-drug regimen for resistant hypertension.
    """
    # Based on the analysis, the following three medications are chosen.
    medication_1 = "Amlodipine"
    medication_2 = "Olmesartan"
    medication_3 = "Chlorthalidone"

    print("For patient MM, the following 3-drug regimen is recommended to maximize hypertension treatment:")
    
    # The final output is structured as a numbered list.
    print(f"1. {medication_1}")
    print(f"2. {medication_2}")
    print(f"3. {medication_3}")
    
    print("\nRationale for this combination:")
    print(f"- {medication_1}: A calcium channel blocker, which is a highly effective first-line agent for African American patients.")
    print(f"- {medication_2}: An angiotensin II receptor blocker (ARB), which is critical for kidney protection in patients with diabetes and works effectively in combination therapy.")
    print(f"- {medication_3}: A potent thiazide-like diuretic. This specific diuretic is not on the patient's list of medications she is unable to take and is a cornerstone for treating resistant hypertension.")

recommend_htn_medications()