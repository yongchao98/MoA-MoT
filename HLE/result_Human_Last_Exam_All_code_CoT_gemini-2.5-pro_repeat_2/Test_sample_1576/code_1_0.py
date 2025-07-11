def find_best_treatment_for_agitation():
    """
    This function analyzes a clinical scenario of a violent, agitated patient
    and determines the best next step based on established medical guidelines.
    """

    # --- Scenario Analysis ---
    print("Analyzing the clinical situation:")
    print("1. The patient is acutely violent, posing an immediate safety risk.")
    print("2. The initial treatment (5mg IM olanzapine) has failed, indicating severe, refractory agitation.")
    print("3. The primary goal is rapid, safe control of agitation to protect the patient and staff.\n")

    # --- Evaluating the Options ---
    print("Evaluating the proposed interventions:")
    print("- Verbal de-escalation alone is insufficient due to existing physical violence.")
    print("- IV administration is high-risk and impractical in a combative patient.")
    print("- For severe agitation unresponsive to a single medication, combination therapy (an antipsychotic + a benzodiazepine) is the standard of care and more effective than escalating a single drug.\n")

    # --- Conclusion ---
    print("Conclusion:")
    print("The most appropriate action is to administer a robust combination of medications.")
    print("Option E provides a potent and effective combination of 10mg IM olanzapine and 2mg IM lorazepam, which is indicated for this level of severe, treatment-resistant agitation.\n")

    # --- Final Answer ---
    best_option_choice = "E"
    best_option_description = "10mg IM olanzapine + 2mg IM lorazepam"

    print(f"The best next step is: {best_option_choice}. {best_option_description}")

    print("\nThe specific dosages in the chosen treatment plan are:")
    print(10)
    print(2)

# Execute the function to get the answer
find_best_treatment_for_agitation()