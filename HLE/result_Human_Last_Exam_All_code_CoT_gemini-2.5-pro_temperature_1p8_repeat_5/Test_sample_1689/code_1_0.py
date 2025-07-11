def solve_dermatology_case():
    """
    This script analyzes the clinical case to determine the best next diagnostic step.
    """

    # --- Case Data from the text ---
    patient_age = 43
    blood_pressure_systolic = 142
    blood_pressure_diastolic = 83
    heart_rate = 87
    respiratory_rate = 15
    
    # --- Key Findings ---
    primary_symptom = "Itchy rash in both axillae"
    key_history = "Started wearing new workout clothes 3 weeks ago"
    key_physical_finding = "Rash on posterior axillary folds, sparing the vaults"
    
    # --- Logical Analysis ---
    print("Analyzing the case of a {}-year-old man.".format(patient_age))
    print("Vital signs include BP: {}/{}, HR: {}, RR: {}.".format(
        blood_pressure_systolic, blood_pressure_diastolic, heart_rate, respiratory_rate))
    print("\nStep 1: Identify the most likely diagnosis based on the evidence.")
    print(" - The patient's history points to a new exposure: {}.".format(key_history))
    print(" - The rash location (sparing the vault) is characteristic of textile dermatitis, not deodorant dermatitis.")
    print(" - Conclusion: The most likely diagnosis is Allergic Contact Dermatitis from clothing.")

    print("\nStep 2: Evaluate the diagnostic options based on this conclusion.")
    print(" A. Skin Biopsy: Too invasive for an initial step with a clear suspected diagnosis.")
    print(" B. KOH Preparation: Unlikely to be positive as the presentation is not typical for a fungal infection.")
    print(" C. Topical Steroid: This is a treatment, not a step in diagnosis.")
    print(" D. Patch Test: This is the gold-standard test to confirm allergic contact dermatitis and identify the specific allergen.")

    print("\nStep 3: Determine the final answer.")
    print("The best next step to confirm the diagnosis is the Patch Test.")
    
    final_answer_option = 'D'
    print("\nFinal Answer Choice: {}".format(final_answer_option))

solve_dermatology_case()