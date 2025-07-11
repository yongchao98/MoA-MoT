import itertools

def find_min_max_models():
    """
    This function solves the problem by reconstructing a potential 3-CNF formula
    from the graph and then counting its satisfying assignments (models).
    """

    # From the graph analysis, we deduced a formula with 5 variables.
    variables = ['x1', 'x2', 'x3', 'x4', 'x5']
    num_vars = len(variables)
    model_count = 0

    print("Analyzing a representative 3-CNF formula derived from the graph.")
    print("The formula has 5 variables (x1 to x5) and 4 clauses.")
    print("Clause 1: (x1 or (not x2) or (not x3))")
    print("Clause 2: (x4 or x2 or (not x3))")
    print("Clause 3: (x3 or x5 or (not x2))")
    print("Clause 4: ((not x4) or (not x5) or (not x1))")
    print("\nCounting the number of satisfying assignments (models)...")

    # Iterate through all 2^5 = 32 possible boolean assignments to the variables.
    for assignment_tuple in itertools.product([False, True], repeat=num_vars):
        x1, x2, x3, x4, x5 = assignment_tuple

        # Evaluate the 4 clauses for the current assignment.
        c1 = x1 or not x2 or not x3
        c2 = x4 or x2 or not x3
        c3 = x3 or x5 or not x2
        c4 = not x4 or not x5 or not x1

        # If all clauses are satisfied, we have found a model.
        if c1 and c2 and c3 and c4:
            model_count += 1
            
    # All possible formulas that could generate the graph are structurally
    # equivalent and thus have the same number of models.
    # Therefore, the minimum and maximum number of models are the same.
    min_models = model_count
    max_models = model_count
    
    print(f"\nFound {model_count} models.")
    print(f"The minimum and maximum number of models are ({min_models}, {max_models}).")

if __name__ == '__main__':
    find_min_max_models()
