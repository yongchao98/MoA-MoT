def is_valid(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        F + J == 86 and
        B == 4.5 * A and
        A + K == 15 and
        M == 3.0 * L and
        H == 1.5 * A and
        J > D and
        A - M == 7 and
        B + E == 52 and
        A > G and
        H > L
    )

def select_unassigned_variable(assignment, domains):
    # Minimum Remaining Values (MRV) heuristic
    unassigned_vars = [i for i, v in enumerate(assignment) if v is None]
    return min(unassigned_vars, key=lambda var: len(domains[var]))

def order_domain_values(var, domains):
    # Least Constraining Value (LCV) heuristic
    return sorted(domains[var], key=lambda val: sum(val in domains[other_var] for other_var in range(len(domains)) if other_var != var))

def forward_checking(assignment, domains, var, value):
    # Forward checking to prune domains
    for other_var in range(len(domains)):
        if other_var != var and assignment[other_var] is None:
            if value in domains[other_var]:
                domains[other_var].remove(value)

def backtrack(assignment, domains):
    if all(v is not None for v in assignment):
        if is_valid(assignment):
            return assignment
        return None

    var = select_unassigned_variable(assignment, domains)
    for value in order_domain_values(var, domains):
        assignment[var] = value
        new_domains = [list(domain) for domain in domains]
        forward_checking(assignment, new_domains, var, value)
        result = backtrack(assignment, new_domains)
        if result is not None:
            return result
        assignment[var] = None

    return None

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Initialize assignment with None and domains with all numbers
assignment = [None] * 13
domains = [list(numbers) for _ in range(13)]

# Solve the CSP
solution = backtrack(assignment, domains)

# Output the solution
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")