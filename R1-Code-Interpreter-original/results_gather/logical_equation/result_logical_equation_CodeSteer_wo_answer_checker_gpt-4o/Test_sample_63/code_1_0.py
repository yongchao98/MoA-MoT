from itertools import permutations
from multiprocessing import Pool

# Define the numbers and the letters
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Preprocess constraints to reduce possibilities
def preprocess_constraints():
    possibilities = {letter: set(numbers) for letter in 'ABCDEFGHIJKLM'}
    
    # Apply direct constraints
    possibilities['E'] = {x - 4 for x in possibilities['A']}
    possibilities['J'] = {x / 4 for x in possibilities['I']}
    possibilities['B'] = {x - 14 for x in possibilities['D']}
    possibilities['M'] = {26 - x for x in possibilities['H']}
    possibilities['F'] = {x / 2 for x in possibilities['A']}
    possibilities['H'] = {x + 4 for x in possibilities['I']}
    
    return possibilities

# Check if the current assignment is valid
def is_valid(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        A - E == 4 and
        I == 4 * J and
        D - B == -14 and
        H + M == 26 and
        L == 3.2 * A and
        L == 1.6 * F and
        F + L == 26 and
        H - I == -4 and
        D == 1.5 * H
    )

# Use backtracking with constraint propagation
def backtrack(assignment, possibilities):
    if len(assignment) == 13:
        if is_valid(assignment):
            return assignment
        return None

    # Select the next variable to assign using MRV heuristic
    unassigned_vars = [var for var in 'ABCDEFGHIJKLM' if var not in assignment]
    next_var = min(unassigned_vars, key=lambda var: len(possibilities[var]))

    for value in possibilities[next_var]:
        if value not in assignment:
            new_assignment = assignment + [value]
            result = backtrack(new_assignment, possibilities)
            if result is not None:
                return result
    return None

# Parallelized search function
def parallel_search(possibilities):
    return backtrack([], possibilities)

# Start the backtracking process with parallelization
if __name__ == '__main__':
    possibilities = preprocess_constraints()
    with Pool() as pool:
        results = pool.map(parallel_search, [possibilities])
        for solution in results:
            if solution:
                print(f"<<<{solution}>>>")
                break