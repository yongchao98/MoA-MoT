def solve_medical_case():
    """
    Analyzes the clinical case and determines the best next diagnostic step.
    """
    # Patient data from the case study
    patient_age = 43
    blood_pressure = "142/83 mm Hg"
    heart_rate = 87
    rash_onset_weeks_ago = 1
    new_activity_start_weeks_ago = 3
    key_symptom = "rash that affects both axillae"
    key_finding = "sparing of axillary vaults"

    print("Step 1: Analyzing the patient's history.")
    print(f"A {patient_age}-year-old man presents with a '{key_symptom}'.")
    print(f"The rash started {rash_onset_weeks_ago} week ago, and the patient began a new workout program {new_activity_start_weeks_ago} weeks ago, which likely involves new workout clothing.")
    print("-" * 20)

    print("Step 2: Analyzing the physical examination findings.")
    print(f"The patient's blood pressure is {blood_pressure} and heart rate is {heart_rate} beats/min.")
    print(f"The most critical skin finding is the distribution of the rash, which involves the posterior axillary folds with '{key_finding}'.")
    print("This specific location strongly suggests an irritant or allergen in the patient's clothing (textile dermatitis) rather than a product applied directly to the armpit, like deodorant.")
    print("-" * 20)

    print("Step 3: Evaluating the diagnostic options.")
    print("The suspected diagnosis is Allergic Contact Dermatitis from textiles.")
    print("A. Skin biopsy: Too invasive for the initial workup of this suspected condition.")
    print("B. KOH preparation: Tests for fungus, which is less likely given the clinical picture.")
    print("C. Topical steroid: This is a treatment to reduce inflammation, not a diagnostic test to find the cause.")
    print("D. Patch test: This is the gold-standard test to identify the specific allergen causing allergic contact dermatitis. It directly tests the hypothesis derived from the history and physical exam.")
    print("-" * 20)
    
    print("Step 4: Conclusion.")
    print("The best next step to confirm the diagnosis and identify the specific trigger is the Patch test.")

    # The final answer in the specified format
    print("<<<D>>>")

solve_medical_case()