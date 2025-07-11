def find_best_next_step():
    """
    Analyzes a clinical scenario to determine the best next step for an agitated patient.
    """

    # 1. Define the patient's current state
    patient_is_violent = True
    initial_drug = "olanzapine"
    initial_dose_mg = 5
    initial_treatment_failed = True

    print("Step 1: Assessing the current situation.")
    if patient_is_violent:
        print("Patient is actively violent. The priority is immediate safety for the patient and staff.")
    print(f"An initial dose of {initial_dose_mg}mg IM {initial_drug} was given with no improvement.")
    print("-" * 30)

    # 2. Evaluate verbal de-escalation (Option C)
    print("Step 2: Evaluating verbal de-escalation.")
    if patient_is_violent:
        print("Option C (Verbal de-escalation) is inappropriate as a *sole next step* due to ongoing violence. Safety must be established first.")
    else:
        print("Verbal de-escalation would be a primary option.")
    print("-" * 30)

    # 3. Evaluate pharmacological options based on route and efficacy
    print("Step 3: Evaluating pharmacological options.")
    print("Option A (2mg IV lorazepam) is high-risk. Obtaining IV access on a violent patient is dangerous.")
    print("\nConsidering the remaining IM options:")
    if initial_treatment_failed:
        print("The initial monotherapy (single drug) failed. A more potent intervention is required.")
        print("Option E combines a higher dose of the antipsychotic with a benzodiazepine, which is a standard, synergistic, and highly effective approach for severe agitation.")

    # 4. Final Recommendation
    final_olanzapine_dose = 10
    final_lorazepam_dose = 2
    best_option = "E"
    
    print("-" * 30)
    print("Conclusion: The most appropriate next step is a combination of medications for rapid and effective control.")
    print(f"The recommended choice is Option {best_option}: 10mg IM olanzapine + 2mg IM lorazepam.")
    
    print("\nFinal Recommended Doses:")
    print(f"Olanzapine dose = {final_olanzapine_dose} mg")
    print(f"Lorazepam dose = {final_lorazepam_dose} mg")

# Run the analysis
find_best_next_step()