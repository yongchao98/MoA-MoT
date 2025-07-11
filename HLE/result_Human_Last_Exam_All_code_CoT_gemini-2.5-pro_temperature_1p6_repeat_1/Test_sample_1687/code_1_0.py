def diagnose_case():
    """
    This script analyzes the clinical vignette to determine the most likely diagnosis
    by breaking down the key findings and ruling out alternatives.
    """
    # Key numerical data from the case study
    initial_hb = 11.7
    second_hb = 6.5
    final_heart_rate = 105
    final_bp_systolic = 100
    final_bp_diastolic = 65

    print("Step 1: Analyze the patient's history and procedure.")
    print("The patient had a 'slightly difficult' colonoscopy, which is a known risk factor for iatrogenic splenic injury due to traction on the splenocolic ligament. Critically, the report states 'No polypectomy was performed,' which rules out postpolypectomy syndrome.")
    print("-" * 30)

    print("Step 2: Evaluate the classic signs and symptoms.")
    print("The patient presented with a constellation of symptoms highly suggestive of splenic injury:")
    print(" - Location of Pain: 'Upper abdominal pain' and 'tenderness to palpation involving her left upper quadrant'.")
    print(" - Referred Pain: 'Left-shoulder discomfort', a classic finding known as Kehr's sign, caused by diaphragmatic irritation from blood.")
    print("-" * 30)

    print("Step 3: Assess evidence of hemorrhage.")
    print("The patient's clinical deterioration points directly to massive, acute bleeding:")
    print(f" - The hemoglobin level dropped precipitously from {initial_hb} g/dL to {second_hb} g/dL.")
    print(" - This blood loss led to hypovolemic shock, evidenced by the development of tachycardia and hypotension.")
    print("-" * 30)
    
    print("Step 4: Formulate a final diagnostic 'equation' based on the key findings.")
    print("The final diagnosis can be summarized with the following clinical markers:")
    print(f"(Difficult Colonoscopy) + (LUQ Pain) + (Left Shoulder Pain) + (Hemoglobin Drop: {initial_hb} -> {second_hb}) + (Shock: HR {final_heart_rate}, BP {final_bp_systolic}/{final_bp_diastolic})")
    print("\nThis combination of findings makes Splenic Laceration the most likely diagnosis.")

diagnose_case()
<<<C>>>