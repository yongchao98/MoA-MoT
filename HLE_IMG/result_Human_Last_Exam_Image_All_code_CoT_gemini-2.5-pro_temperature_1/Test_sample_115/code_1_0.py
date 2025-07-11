from itertools import product

def count_models(num_vars, clauses):
    """Counts the number of satisfying assignments for a CNF formula."""
    model_count = 0
    # Generate all 2^n possible truth assignments for the variables
    assignments = product([False, True], repeat=num_vars)

    for assignment in assignments:
        # A dictionary to map variable index to its truth value
        var_map = {i + 1: val for i, val in enumerate(assignment)}
        
        is_satisfied = True
        # Check if the current assignment satisfies all clauses
        for clause in clauses:
            clause_val = False
            for literal in clause:
                var_index = abs(literal)
                is_negated = literal < 0
                
                # XOR for negation
                if var_map[var_index] ^ is_negated:
                    clause_val = True
                    break
            
            if not clause_val:
                is_satisfied = False
                break
        
        if is_satisfied:
            model_count += 1
            
    return model_count

def solve():
    """
    Calculates the min and max number of models for the possible formulas
    derived from the graph structure.
    """
    
    # Base clauses from the 3-triangle substructure.
    # Variables are x1, x2, x3, x4. Literals are represented by signed integers.
    # C1 = (x4 | x1 | x3)
    # C3 = (~x2 | ~x4 | ~x3)
    c1 = {4, 1, 3}
    c3 = {-2, -4, -3}
    
    # Case 1: Max models. l_4 is a new variable x5.
    # C2 = (~x1 | x5 | x2)
    clauses_max = [c1, {-1, 5, 2}, c3]
    num_vars_max = 5
    max_models = count_models(num_vars_max, clauses_max)

    # Case 2: Min models exploration
    # Subcase 2a: l_4 makes a clause tautological. Let l_4 = x1.
    # The second clause (~x1 | x1 | x2) is always true and can be removed.
    clauses_min_a = [c1, c3]
    num_vars_min_a = 4 # Only variables x1, x2, x3, x4 matter
    min_models_a = count_models(num_vars_min_a, clauses_min_a)

    # Subcase 2b: l_4 is an existing variable. Let l_4 = x3.
    # C2 = (~x1 | x3 | x2)
    clauses_min_b = [c1, {-1, 3, 2}, c3]
    num_vars_min_b = 4
    min_models_b = count_models(num_vars_min_b, clauses_min_b)

    # The minimum is the smaller of the two cases explored
    min_models = min(min_models_a, min_models_b)
    
    print(f"Maximum number of models: {max_models}")
    print(f"Minimum number of models: {min_models}")
    print(f"The answer is the pair (min, max).")

solve()

# The final answer format
print("<<<(10, 20)>>>")