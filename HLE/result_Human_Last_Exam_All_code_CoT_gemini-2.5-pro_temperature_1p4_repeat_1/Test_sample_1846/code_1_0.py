def analyze_tnbc_treatment_efficacy():
    """
    Analyzes and prints the findings from clinical trials on PD-1 inhibitors
    for Triple Negative Breast Cancer (TNBC).
    """

    # Landmark clinical trials, such as KEYNOTE-355 (testing the PD-1 inhibitor pembrolizumab),
    # have evaluated the addition of immunotherapy to chemotherapy in TNBC.
    # The effectiveness was measured in different patient groups.

    # Group 1: Intention-to-treat (ITT) population. This includes all randomized patients,
    # regardless of their PD-L1 expression status.
    itt_population_result = "In the ITT population, the improvement in overall survival did not reach statistical significance."

    # Group 2: PD-L1-positive population. This subgroup includes patients whose tumors
    # express the PD-L1 biomarker, often measured by a Combined Positive Score (CPS).
    pdl1_positive_result = "In the PD-L1-positive population (specifically those with a higher CPS score), a statistically significant and clinically meaningful prolongation of overall survival was observed."
    
    # Group 3: PD-L1-negative population. This subgroup includes patients whose tumors do not
    # express the PD-L1 biomarker.
    pdl1_negative_result = "In the PD-L1-negative population, there was no significant benefit in overall survival."

    print("Analysis of PD-1 Inhibitor Efficacy in Triple Negative Breast Cancer:")
    print("="*70)

    print("\nPopulation Group: Intention-to-treat (All Patients)")
    print(f"Finding: {itt_population_result}")

    print("\nPopulation Group: PD-L1-positive")
    print(f"Finding: {pdl1_positive_result}")

    print("\nPopulation Group: PD-L1-negative")
    print(f"Finding: {pdl1_negative_result}")

    print("\n" + "="*70)
    print("Conclusion: The addition of PD-1 inhibitors to chemotherapy demonstrates a prolonged overall survival benefit specifically in the PD-L1-positive population, not in the overall intention-to-treat or PD-L1-negative populations.")

# Execute the analysis
analyze_tnbc_treatment_efficacy()