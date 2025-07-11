def analyze_population_genetics_scenario():
    """
    Analyzes the provided population genetics scenario to determine which statements
    must always be true.
    """

    print("Analyzing the truth of four statements based on the given population parameters...")
    print("-" * 70)

    # --- Conditions Summary ---
    # H_squared = 0.7, meaning Phenotypic Variance Vp = Vg + Ve, and Ve is 30% of Vp.
    # No selection, no mutation, no drift (infinite pop), no migration.
    # Random mating throughout the entire region.

    # --- Statement 1 Analysis ---
    statement_1_true = True
    print("Statement 1: 'There is no selection occurring on the phenotype measured.'")
    print("Analysis: The prompt explicitly states the trait has 'no bearing on fitness' and 'all genotypes have equal fitness'. This is the definition of no selection.")
    print("Conclusion: This statement MUST be true.\n")

    # --- Statement 2 Analysis ---
    statement_2_true = False
    print("Statement 2: 'Parents will not raise their offspring.'")
    print("Analysis: The prompt gives no information on parental behavior. Parental care is an environmental factor. Since broad-sense heritability is 0.7, environmental variance exists and could include effects of parental care. This statement is not required by the given information.")
    print("Conclusion: This statement does NOT have to be true.\n")

    # --- Statement 3 Analysis ---
    statement_3_true = True
    print("Statement 3: 'The population will never speciate even in future generations.'")
    print("Analysis: Speciation requires evolution. The prompt has systematically removed all mechanisms of evolution: no selection, no mutation, no genetic drift, and no barriers to gene flow (due to random mating). Within the defined constraints of the problem, the population cannot evolve or speciate.")
    print("Conclusion: This statement MUST be true.\n")

    # --- Statement 4 Analysis ---
    statement_4_true = False
    print("Statement 4: 'The researcher will not find a substantial difference in the phenotype measured between the west and east groups of individuals.'")
    print("Analysis: While random mating ensures no genetic differences between east and west, the phenotype is also 30% determined by the environment (since H^2 = 0.7). The problem does not state the environment is uniform. A systematic environmental difference (e.g., climate, diet) between east and west could cause a substantial difference in the measured phenotype.")
    print("Conclusion: This statement does NOT have to be true.\n")

    # --- Final Conclusion ---
    print("-" * 70)
    print("Summary of conclusions:")
    print("Statement 1 must be true.")
    print("Statement 3 must be true.")
    print("\nThe correct option combines statements 1 and 3 only.")
    print("\nFinal Answer Choice: F")


# Execute the analysis function
analyze_population_genetics_scenario()

print("<<<F>>>")