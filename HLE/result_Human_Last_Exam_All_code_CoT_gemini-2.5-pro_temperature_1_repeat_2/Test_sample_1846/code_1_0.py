def analyze_tnbc_treatment_outcomes():
    """
    Analyzes clinical trial data to determine which population group benefits
    most from PD-1 inhibitors in Triple Negative Breast Cancer (TNBC).
    """

    # Clinical trial findings for overall survival (OS) when adding PD-1 inhibitors to chemotherapy
    outcomes = {
        "Intention-to-treat population": {
            "description": "This group includes all patients in the trial.",
            "os_benefit": "Not consistently significant across major trials for metastatic disease."
        },
        "PD-L1-positive population": {
            "description": "Patients whose tumors express the PD-L1 biomarker.",
            "os_benefit": "Statistically significant and clinically meaningful prolonged overall survival."
        },
        "PD-L1-negative population": {
            "description": "Patients whose tumors do not express the PD-L1 biomarker.",
            "os_benefit": "No significant overall survival benefit observed."
        }
    }

    print("Analyzing Overall Survival Benefit of PD-1 Inhibitors in TNBC by Population Group:")
    print("="*80)

    best_population = None
    max_benefit_description = ""

    for population, data in outcomes.items():
        print(f"Population Group: {population}")
        print(f"Description: {data['description']}")
        print(f"Overall Survival (OS) Finding: {data['os_benefit']}")
        print("-" * 80)
        if "Statistically significant" in data['os_benefit']:
            best_population = population
            max_benefit_description = data['os_benefit']

    print("\nConclusion:")
    print(f"Based on pivotal clinical trials like KEYNOTE-355, the group that demonstrated a {max_benefit_description.lower()}")
    print(f"is the '{best_population}'.")


analyze_tnbc_treatment_outcomes()