def analyze_tnbc_treatment_survival():
    """
    Analyzes and prints the findings from clinical trials regarding PD-1 inhibitors
    in Triple Negative Breast Cancer (TNBC) to determine which population
    group shows prolonged overall survival.
    """

    # Define the patient populations being compared
    population_a = "Intention-to-treat population"
    population_b = "PD-L1-positive population"
    population_d = "PD-L1-negative population"

    print("Analysis of Overall Survival in TNBC Patients Treated with PD-1 Inhibitors + Chemotherapy")
    print("="*80)
    
    # Explain the findings for the PD-L1-positive population
    print(f"Population Group: {population_b}")
    print("Finding: Major clinical trials, such as KEYNOTE-355, have shown that adding a PD-1 inhibitor (like pembrolizumab) to chemotherapy results in a statistically significant and clinically meaningful improvement in overall survival for patients with PD-L1-positive metastatic TNBC.")
    print("This benefit is the basis for the approval and use of this combination in this specific subgroup.")
    print("-" * 80)

    # Explain the findings for the intention-to-treat (ITT) population
    print(f"Population Group: {population_a}")
    print("Finding: In the overall intention-to-treat (ITT) population, which includes both PD-L1-positive and PD-L1-negative patients, a statistically significant improvement in overall survival was NOT observed.")
    print("This indicates the benefit is concentrated in the PD-L1-positive subgroup and is diluted across the entire population.")
    print("-" * 80)
    
    # Explain the findings for the PD-L1-negative population
    print(f"Population Group: {population_d}")
    print("Finding: Patients with PD-L1-negative tumors did not show a significant survival benefit from the addition of a PD-1 inhibitor.")
    print("-" * 80)
    
    # Final conclusion based on the evidence
    print("Conclusion:")
    print(f"The evidence clearly indicates that the treatment with PD-1 inhibitors presents a prolonged overall survival specifically in the {population_b}.")

# Run the analysis
analyze_tnbc_treatment_survival()