from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Constraints
def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I = assignment
    return (
        D + H == 36 and
        G - I == -79 and
        C - B == -20 and
        D == 3 * E and
        B == 4 * D and
        I > C and
        B + C == 52
    )

# Find the correct assignment
for perm in permutations(values):
    if check_constraints(perm):
        print(list(perm))
        break