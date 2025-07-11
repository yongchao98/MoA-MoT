def analyze_statements():
    """
    Analyzes four statements about a hypothetical population based on a given description.
    """

    print("Analyzing which statement must always be true.")
    print("--------------------------------------------------")

    # Statement 1 Analysis
    print("Statement 1: 'There is no selection occurring on the phenotype measured.'")
    print("The problem states 'all genotypes have equal fitness' and the phenotype has 'no bearing on fitness'.")
    print("This is the definition of no selection. So, statement 1 is TRUE.\n")

    # Statement 2 Analysis
    print("Statement 2: 'Parents will not raise their offspring.'")
    print("The problem provides no information about parental care. 'Non-overlapping generations' does not forbid it.")
    print("Therefore, we cannot conclude this must be true. So, statement 2 is FALSE.\n")

    # Statement 3 Analysis
    print("Statement 3: 'The population will never speciate even in future generations.'")
    print("The conditions are described as 'currently true'. They might change in the future.")
    print("If conditions change (e.g., a mountain appears), speciation could occur. 'Never' is too strong a claim.")
    print("Therefore, this does not have to be true. So, statement 3 is FALSE.\n")

    # Statement 4 Analysis
    print("Statement 4: 'The researcher will not find a substantial difference in the phenotype measured between the west and east groups of individuals.'")
    print("Random mating makes the *genetics* uniform. However, heritability is 0.7, meaning environment accounts for 30% of the variance.")
    print("An environmental difference between the west and east could cause a phenotypic difference.")
    print("Therefore, this does not have to be true. So, statement 4 is FALSE.\n")

    print("--------------------------------------------------")
    print("Conclusion: Only statement 1 must always be true.")
    print("This corresponds to option 'A. 1 only'.")

    # The final answer is the choice corresponding to "1 only".
    final_answer = 'A'
    print(f"\nFinal choice: {final_answer}")
    print(f"<<<{final_answer}>>>")

analyze_statements()