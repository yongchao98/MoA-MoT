def find_most_problematic_clause():
    """
    Analyzes terms of service clauses to identify the one most likely to be a
    problematic contract of adhesion or hide a material term.
    """
    # A score is assigned based on the "surprise factor" of the clause.
    # A higher score means the term is less expected by a reasonable user.
    # 1 = Standard/Expected Term
    # 10 = Highly Unusual/Unexpected Term
    clause_scores = {
        'A': 2, # Standard non-compete/benchmarking clause.
        'B': 1, # Standard content license needed for a service to function.
        'C': 1, # Standard late fees and non-refundable terms.
        'D': 1, # Standard anti-scraping clause.
        'E': 10, # Contains a highly unusual and specific geographic restriction (Illinois).
        'F': 3, # A modern, but increasingly standard, anti-AI-training clause.
        'G': 5  # Controversial social advertising clause, but well-established on major platforms.
    }

    # "Equation" Part: Print the evaluation of each term's score.
    # This fulfills the prompt's requirement to show numbers in an equation-like format.
    print("Evaluating clauses based on a 'surprise factor' score (1-10):")
    for clause, score in clause_scores.items():
        print(f"Clause {clause} has a surprise score of: {score}")

    # Find the clause with the maximum score.
    most_problematic_clause = max(clause_scores, key=clause_scores.get)
    max_score = clause_scores[most_problematic_clause]

    print(f"\nThe clause with the highest surprise factor score is '{most_problematic_clause}' with a score of {max_score}.")
    print("\nThis clause is the most likely to be considered a contract of adhesion with hidden material terms because it contains a highly specific and unexpected restriction (forbidding research on individuals in Illinois) that a reasonable user would not anticipate.")
    print("\nFinal Answer:")
    print(most_problematic_clause)

find_most_problematic_clause()