def find_benefiting_population():
    """
    Analyzes clinical trial data to determine which TNBC population shows
    prolonged overall survival with PD-1 inhibitors plus chemotherapy.
    """
    # This data structure represents simplified findings from key clinical trials
    # like KEYNOTE-355.
    # 'Benefit' refers to the comparison of (PD-1 inhibitor + chemo) vs. (chemo alone).
    clinical_trial_outcomes = {
        "A. Intention-to-treat population": "No statistically significant prolonged overall survival benefit.",
        "B. PD-L1-positive population": "Statistically significant prolonged overall survival benefit observed.",
        "D. PD-L1-negative population": "No prolonged overall survival benefit observed."
    }

    benefiting_group = None
    explanation = ""

    print("Analyzing clinical trial outcomes for PD-1 inhibitors in TNBC...")
    print("-" * 60)

    for population, outcome in clinical_trial_outcomes.items():
        if "prolonged overall survival benefit observed" in outcome.lower():
            benefiting_group = population
            explanation = outcome
            break # Found the primary group with the benefit

    if benefiting_group:
        print(f"Population with proven benefit: {benefiting_group}")
        print(f"Finding: {explanation}")
        print("\nThis leads to the following conclusion:")
        # The final "equation" demonstrates the relationship
        print("Final Equation: Overall_Survival(PD-L1-positive + PD-1_inhibitor) > Overall_Survival(PD-L1-positive + Chemotherapy_alone)")
    else:
        print("No population found with a clear overall survival benefit.")

find_benefiting_population()