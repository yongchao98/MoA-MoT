from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Initial domains for each letter
domains = {
    'A': set(numbers),
    'B': set(numbers),
    'C': set(numbers),
    'D': set(numbers),
    'E': set(numbers),
    'F': set(numbers),
    'G': set(numbers),
    'H': set(numbers),
    'I': set(numbers),
    'J': set(numbers),
    'K': set(numbers),
    'L': set(numbers),
    'M': set(numbers)
}

# Constraints as functions
def apply_constraints(domains):
    # Apply constraints to reduce domains
    for num in numbers:
        if num in domains['B']:
            domains['A'] = {4 * num}
        if num in domains['L']:
            domains['D'] = {3 * num}
            domains['J'] = {num // 3}
        if num in domains['F']:
            domains['B'] = {num - 22}
            domains['C'] = {num + 14}
        if num in domains['I']:
            domains['D'] = {num + 17}
            domains['E'] = {num - 21}
            domains['J'] = {33 - num}
        if num in domains['D']:
            domains['C'] = {num - 29}
            domains['L'] = {num // 3}
        if num in domains['F'] and num in domains['L']:
            if num + domains['L'].pop() == 17:
                domains['L'] = {17 - num}

# Backtracking function
def backtrack(assignment, domains):
    if len(assignment) == 13:
        return assignment

    # Select the next variable to assign
    unassigned_vars = [var for var in domains if var not in assignment]
    var = min(unassigned_vars, key=lambda v: len(domains[v]))

    for value in domains[var]:
        # Create a new assignment with the current value
        new_assignment = assignment.copy()
        new_assignment[var] = value

        # Create new domains with the current value removed
        new_domains = {k: (v - {value} if k != var else {value}) for k, v in domains.items()}

        # Apply constraints to the new domains
        apply_constraints(new_domains)

        # Recursively backtrack with the new assignment and domains
        result = backtrack(new_assignment, new_domains)
        if result is not None:
            return result

    return None

# Find the solution using backtracking with constraint propagation
solution = backtrack({}, domains)

# Print the solution
if solution:
    # Order the solution alphabetically
    ordered_solution = [solution[letter] for letter in sorted(solution.keys())]
    print(f"<<<{ordered_solution}>>>")
else:
    print("No solution found.")