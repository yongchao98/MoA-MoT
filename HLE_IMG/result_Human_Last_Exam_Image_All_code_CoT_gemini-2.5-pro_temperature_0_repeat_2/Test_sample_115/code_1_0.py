import itertools

def evaluate_literal(lit, assignment):
    """
    Evaluates a literal's truth value under a given assignment.
    A literal is a tuple (variable_index, required_polarity).
    An assignment is a tuple of booleans representing the truth values of variables.
    """
    var_index, polarity = lit
    return assignment[var_index] == polarity

def evaluate_clause(clause, assignment):
    """Evaluates a clause's truth value."""
    for lit in clause:
        if evaluate_literal(lit, assignment):
            return True
    return False

def evaluate_formula(formula, assignment):
    """Evaluates a formula's truth value."""
    for clause in formula:
        if not evaluate_clause(clause, assignment):
            return False
    return True

def count_models(formula, num_vars):
    """Counts the number of satisfying assignments (models) for a formula."""
    count = 0
    # Iterate through all 2^num_vars possible assignments
    for assignment_tuple in itertools.product([False, True], repeat=num_vars):
        if evaluate_formula(formula, assignment_tuple):
            count += 1
    return count

def solve():
    """
    Reverses the 3-SAT to IS reduction for the given graph and finds the
    min and max number of models for possible source formulas.
    """
    # There are 8 variables in total, indexed 0 through 7.
    # x1-x4 (indices 0-3) are the variables with negation links.
    # x5-x8 (indices 4-7) are the variables for pure literals.
    num_vars = 8

    # Helper functions to create literal representations
    x1_lit = lambda pol: (0, pol)
    x2_lit = lambda pol: (1, pol)
    x3_lit = lambda pol: (2, pol)
    x4_lit = lambda pol: (3, pol)
    
    # Pure literals can be fixed to positive polarity without loss of generality,
    # as flipping them just renames the variable.
    x5_lit = (4, True)
    x6_lit = (5, True)
    x7_lit = (6, True)
    x8_lit = (7, True)

    min_models = float('inf')
    max_models = float('-inf')

    # Iterate through the 2^4 = 16 polarity choices for x1, x2, x3, x4.
    # p_i = True means the base literal is positive (e.g., l1 = x1).
    # p_i = False means the base literal is negative (e.g., l1 = ~x1).
    for p1, p2, p3, p4 in itertools.product([False, True], repeat=4):
        
        # Define the literals for this specific formula based on polarity choices
        l1 = x1_lit(p1)
        l2 = x2_lit(p2)
        l3 = x3_lit(p3)
        l4 = x4_lit(p4)
        
        # Define their negations
        l1_neg = x1_lit(not p1)
        l2_neg = x2_lit(not p2)
        l3_neg = x3_lit(not p3)
        l4_neg = x4_lit(not p4)
        
        # Based on the graph structure, the four clauses are:
        # C1 = (x5 or l1 or l2)
        # C2 = (x6 or x7 or l3)
        # C3 = (x8 or ~l1 or l4)
        # C4 = (~l2 or ~l4 or ~l3)
        formula = [
            [x5_lit, l1, l2],
            [x6_lit, x7_lit, l3],
            [x8_lit, l1_neg, l4],
            [l2_neg, l4_neg, l3_neg]
        ]
        
        models = count_models(formula, num_vars)
        
        if models < min_models:
            min_models = models
        if models > max_models:
            max_models = models

    print(f"({min_models}, {max_models})")
    # The final answer format required by the prompt
    print(f"<<<({min_models}, {max_models})>>>")

solve()