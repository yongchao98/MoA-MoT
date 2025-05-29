from itertools import permutations

# Define the numbers and the constraints
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

def satisfies_constraints(assignment):
    A, B, C, D, E, F, G, H, I = assignment
    return (
        D == 3 * B and
        A + F == 52 and
        D - H == -71 and
        B + I == 30 and
        D + E == 13 and
        I == 3 * D and
        F == 4 * E
    )

# Try all permutations of the numbers
for perm in permutations(numbers):
    if satisfies_constraints(perm):
        print(f"<<<{list(perm)}>>>")
        break