def solve_population_genetics_problem():
    """
    Analyzes the given statements about a hypothetical population
    to determine which must be true.
    """

    # --- Analysis of each statement ---

    # Statement 1: "There is no selection occurring on the phenotype measured."
    # The prompt explicitly states the phenotype has "no bearing on fitness" and "all genotypes have equal fitness".
    # This is the definition of a scenario with no natural selection.
    analysis_1 = "Statement 1 is TRUE. The problem directly states that the phenotype and associated genotypes have no effect on fitness, which means there is no selection."

    # Statement 2: "Parents will not raise their offspring."
    # The prompt mentions "discrete and non-overlapping generations". This is a life-cycle model and does not
    # forbid parental care, which would be an environmental factor. Since H^2 = 0.7, there is a 30% environmental
    # component to phenotypic variance, which could include parental care.
    analysis_2 = "Statement 2 is NOT NECESSARILY TRUE. 'Non-overlapping generations' is a modeling assumption about reproductive timing and does not rule out parental care as an environmental influence."

    # Statement 3: "The population will never speciate even in future generations."
    # The conditions described are "currently true". They do not apply indefinitely into the future.
    # Future events like mutations, environmental changes leading to selection, or the rise of a geographic
    # barrier could lead to speciation. The word "never" is too absolute.
    analysis_3 = "Statement 3 is NOT NECESSARILY TRUE. The conditions preventing evolution are only for the present. They could change in the future, allowing speciation to occur."

    # Statement 4: "The researcher will not find a substantial difference in the phenotype measured between the west and east groups of individuals."
    # The population has "random mating", which implies panmixia. This means individuals mate regardless of their
    # geographic location, leading to a single, homogeneous gene pool across the entire region.
    # With no genetic or systematic environmental differences between the west and east, their average phenotypes will be the same.
    analysis_4 = "Statement 4 is TRUE. Random mating (panmixia) across the whole population prevents genetic differentiation between geographic sub-regions. Therefore, the average phenotype in the west and east will be identical."

    # --- Conclusion ---
    conclusion = "Based on the analysis, only statements 1 and 4 must always be true. This corresponds to option G."

    print("Step-by-step analysis:")
    print(f"1. {analysis_1}")
    print(f"2. {analysis_2}")
    print(f"3. {analysis_3}")
    print(f"4. {analysis_4}")
    print("\n" + conclusion)

solve_population_genetics_problem()