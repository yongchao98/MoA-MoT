def analyze_agitation_treatment():
    """
    Analyzes the clinical scenario of an agitated patient and determines the best next step.
    """
    initial_olanzapine_mg = 5
    patient_status = "Violent, no improvement after initial medication."

    print("Clinical Scenario Analysis:")
    print(f"Initial State: A patient is acutely agitated and has already received {initial_olanzapine_mg} mg of olanzapine (Zyprexa) IM with no improvement.")
    print("The primary goal is to gain rapid, safe control of the agitation to ensure patient and staff safety.")
    print("-" * 30)
    print("Evaluation of Next Steps:")
    print("1. Verbal de-escalation alone is insufficient due to active violence.")
    print("2. IV administration is impractical and dangerous in a combative patient.")
    print("3. Combination therapy (antipsychotic + benzodiazepine) is often more effective than increasing the dose of a single agent that has already failed.")
    print("4. High-dose combinations increase risks of adverse effects like respiratory depression, especially with unknown patient history.")
    print("-" * 30)

    # The best choice is B: 2mg IM lorazepam + 5mg olanzapine IM
    best_choice = "B. 2mg IM lorazepam + 5mg olanzapine IM"
    added_lorazepam_mg = 2
    added_olanzapine_mg = 5

    total_olanzapine_mg = initial_olanzapine_mg + added_olanzapine_mg
    total_lorazepam_mg = added_lorazepam_mg

    print(f"Recommended Next Step: {best_choice}")
    print("This approach adds a synergistic agent (lorazepam) and increases the olanzapine to a standard, effective dose without being excessive.")
    print("\nFinal Dosage Calculation:")
    print(f"Initial Olanzapine Dose: {initial_olanzapine_mg} mg")
    print(f"Added Medication: {added_lorazepam_mg} mg Lorazepam + {added_olanzapine_mg} mg Olanzapine")
    print("=" * 25)
    print(f"Total Administered Dose: {total_olanzapine_mg} mg Olanzapine + {total_lorazepam_mg} mg Lorazepam")

analyze_agitation_treatment()