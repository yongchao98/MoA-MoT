from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Constraints as functions
def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    return (
        D == 3.0 * I and
        C + H == 20 and
        A > F and
        F + G == 38 and
        L - G == -26 and
        C + F == 17 and
        E + L == 60 and
        M == 4.5 * L and
        A - K == -17 and
        D + H == 8
    )

# Iterate over all permutations and check constraints
for perm in permutations(numbers):
    if check_constraints(perm):
        print(f"<<<{list(perm)}>>>")
        break