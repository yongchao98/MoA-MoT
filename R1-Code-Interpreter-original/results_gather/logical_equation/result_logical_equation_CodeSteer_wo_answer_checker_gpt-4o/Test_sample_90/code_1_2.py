def solve_csp():
    from itertools import permutations

    # Given numbers
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

    # Initial domains for each variable
    domains = {var: set(numbers) for var in "ABCDEFGHIJKLM"}

    # Constraints
    constraints = [
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: E + K == 81,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: K > J,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: E == 1.5 * G,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: B + J == 5,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: M - F == -6,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: K == 4.5 * M,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: M > B,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: F == 3.2 * C,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: C - E == -31,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: L == 1.4 * C,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: A + M == 38,
        lambda A, B, C, D, E, F, G, H, I, J, K, L, M: F + I == 112,
    ]

    # Backtracking with forward checking
    def backtrack(assignment):
        if len(assignment) == 13:
            A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
            if all(constraint(A, B, C, D, E, F, G, H, I, J, K, L, M) for constraint in constraints):
                print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}, {L}, {M}]>>>")
                return True
            return False

        # Select the next variable to assign
        unassigned_vars = [var for var in "ABCDEFGHIJKLM" if var not in assignment]
        var = unassigned_vars[0]

        # Try each value in the domain of the variable
        for value in domains[var]:
            if value not in assignment:
                new_assignment = assignment + [value]
                if backtrack(new_assignment):
                    return True

        return False

    # Start backtracking with an empty assignment
    backtrack([])

solve_csp()