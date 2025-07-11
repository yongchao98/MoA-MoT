def solve_geology_task():
    """
    Analyzes geological statements about the North American Cordillera,
    classifies them as Consensus (C) or Debated (D), and prints the result.
    """
    
    # Each tuple contains the statement number and its classification.
    # C = Consensus, D = Debated
    classifications = [
        (1, "C"),
        (2, "D"),
        (3, "C"),
        (4, "D"),
        (5, "C"),
        (6, "C"),
        (7, "C"),
        (8, "D"),
        (9, "D"),
        (10, "D")
    ]
    
    print("Classification for each numbered statement:")
    
    final_string = ""
    for num, classification in classifications:
        print(f"Statement {num}: {classification}")
        final_string += classification
        
    print("\n" + "="*30)
    print("Final Continuous String Answer:")
    print(final_string)
    print("="*30)

solve_geology_task()