import itertools

def count_models_for_formula(clauses):
    """Counts models for a CNF formula represented as a list of clauses."""
    if not clauses:
        # Should determine the number of variables if there's no clause, 
        # but in this problem that won't happen.
        # This is a bit ambiguous, but let's assume vars are inferred.
        return 0

    variables = set()
    for clause in clauses:
        for literal in clause:
            variables.add(abs(literal))
    
    if not variables:
        # If there are clauses but no variables, e.g. [[]], it's unsat.
        # If clauses is empty, it's a tautology, but what's the var set?
        # Let's assume this path isn't taken.
        return 0

    var_list = sorted(list(variables))
    var_map = {var: i for i, var in enumerate(var_list)}
    num_vars = len(var_list)
    count = 0

    for values in itertools.product([False, True], repeat=num_vars):
        assignment = {var: values[var_map[var]] for var in var_list}
        
        formula_is_true = True
        for clause in clauses:
            clause_is_true = False
            for literal in clause:
                var_index = abs(literal)
                required_value = literal > 0
                if assignment[var_index] == required_value:
                    clause_is_true = True
                    break
            if not clause_is_true:
                formula_is_true = False
                break
        
        if formula_is_true:
            count += 1
            
    return count

def solve():
    """
    Analyzes the possible formulas and finds the min and max number of models.
    """
    
    # Base clauses from graph analysis (see plan)
    # C2 = ~x1 V x3 V x4
    # C4 = ~x2 V ~x4 V ~x5
    c2 = [-1, 3, 4]
    c4 = [-2, -4, -5]

    # Different possibilities for the formula based on choices for C1 and C3
    model_counts = []

    # Case 1: L1=x6, L2=x7 (Formula with 7 variables)
    # C1 = x1 V x2 V x6
    # C3 = ~x3 V x5 V x7
    formula_7_vars = [[1, 2, 6], c2, [-3, 5, 7], c4]
    count1 = count_models_for_formula(formula_7_vars)
    model_counts.append(count1)

    # Case 2: One clause is a tautology, the other introduces a new variable (6 vars)
    # 2a: C1 is tautology, C3 introduces x7
    formula_6_vars_a = [c2, [-3, 5, 7], c4]
    count2a = count_models_for_formula(formula_6_vars_a)
    model_counts.append(count2a)

    # 2b: C1 introduces x6, C3 is tautology
    formula_6_vars_b = [[1, 2, 6], c2, c4]
    count2b = count_models_for_formula(formula_6_vars_b)
    model_counts.append(count2b)

    # Case 3: Both C1 and C3 are tautologies (5 vars)
    formula_5_vars = [c2, c4]
    count3 = count_models_for_formula(formula_5_vars)
    model_counts.append(count3)

    min_models = min(model_counts)
    max_models = max(model_counts)

    print(f"The number of models for the 7-variable formula (most constrained) is {count1}.")
    print(f"The number of models for the 6-variable formula (C1 tautology) is {count2a}.")
    print(f"The number of models for the 6-variable formula (C3 tautology) is {count2b}.")
    print(f"The number of models for the 5-variable formula (least constrained) is {count3}.")
    print(f"\nThe set of possible model counts is {sorted(list(set(model_counts)))}.")
    print(f"\nMinimum number of models: {min_models}")
    print(f"Maximum number of models: {max_models}")
    print(f"\nThe pair (min, max) is ({min_models}, {max_models}).")


solve()
<<<24, 76>>>