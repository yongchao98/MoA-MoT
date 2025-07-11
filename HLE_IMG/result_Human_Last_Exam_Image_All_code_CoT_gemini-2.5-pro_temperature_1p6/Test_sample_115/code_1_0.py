import itertools

def count_models(variables, formula_func):
    """
    Counts the number of satisfying assignments (models) for a given boolean formula.

    Args:
        variables (list): A list of variable names, e.g., ['x', 'y', 'z'].
        formula_func (function): A function that takes a dictionary representing an
                                 assignment (e.g., {'x': True, 'y': False}) and
                                 returns the boolean value of the formula.
    Returns:
        int: The number of models.
    """
    num_vars = len(variables)
    model_count = 0
    
    # Generate all possible assignments of True/False to the variables
    # 2**num_vars assignments in total.
    for p in itertools.product([True, False], repeat=num_vars):
        assignment = dict(zip(variables, p))
        if formula_func(assignment):
            model_count += 1
            # print(f"Found model: {assignment}") # Uncomment to see models
            
    return model_count

# Case 1: Corresponds to the formula with the maximum number of models.
# Let the formula be equivalent to MAJ(x,y,z) = (x and y) or (y and z) or (z and x).
# This represents a logical structure derivable from one of the graph's factorizations.
def formula_max(a):
    x, y, z = a['x'], a['y'], a['z']
    return (x and y) or (y and z) or (z and x)

# Case 2: Corresponds to the formula with the minimum number of models.
# Let the formula be equivalent to just 'a'.
# This represents a logical structure from the other possible factorization.
def formula_min(a):
    return a['a']

# --- Calculations ---

# Max models from a formula with 3 variables
max_variables = ['x', 'y', 'z']
max_models = count_models(max_variables, formula_max)

# Min models from a formula with 2 variables
min_variables = ['a', 'b']
min_models = count_models(min_variables, formula_min)

print(f"The minimum number of models is: {min_models}")
print(f"The maximum number of models is: {max_models}")
print(f"The final answer as a pair (min, max) is: ({min_models}, {max_models})")