import sys

def run_diagnostic_analysis():
    """
    Analyzes clinical findings to determine the most likely diagnosis.
    """
    # Patient Data from the case
    age = 64
    bmi = 39
    smoking_history_years = 15

    # Key Clinical Findings
    lesion_locations = ["axillary folds", "inframammary folds", "inguinal regions"]
    lesion_types = ["large bullae", "erythematous plaques", "purulent nodules"]
    
    # --- Diagnostic Reasoning Simulation ---
    print("Step 1: Analyzing Patient Risk Factors")
    # This section fulfills the instruction to output each number from the case.
    print(f" - Patient Age: {age} years")
    print(f" - Body Mass Index (BMI): {bmi} (Obesity is a major risk factor for HS)")
    print(f" - Smoking History: {smoking_history_years} years (Smoking is a major risk factor for HS)")
    
    print("\nStep 2: Evaluating Clinical Presentation")
    print(f" - Lesion Locations: {', '.join(lesion_locations)}")
    print("   -> These are classic intertriginous sites for Hidradenitis Suppurativa (HS).")
    
    print(f" - Lesion Types: {', '.join(lesion_types)}")
    print("   -> The presence of 'purulent nodules' is a hallmark feature of HS.")

    print("\nStep 3: Conclusion")
    print("The combination of significant risk factors (obesity, smoking) with lesions in classic locations")
    print("and the presence of cardinal signs like purulent nodules strongly indicates a diagnosis of")
    print("Hidradenitis Suppurativa over other differentials like psoriasis or dermatitis.")
    
    # Final Answer
    final_diagnosis = "C. Hidradenitis Supportiva"
    print(f"\nMost Likely Diagnosis: {final_diagnosis}")

# Execute the analysis
run_diagnostic_analysis()