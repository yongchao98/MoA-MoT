from itertools import permutations

# Given numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Constraints
def satisfies_constraints(assignment):
    A, B, C, D, E, F, G, H, I = assignment
    return (
        A + F == 141 and
        A + G == 124 and
        I - C == 26 and
        H - G == 33 and
        C == 2 * G and
        D - H == -27 and
        D - G == 6 and
        H - E == -84
    )

# Try all permutations of the numbers
for perm in permutations(numbers):
    if satisfies_constraints(perm):
        # If a valid assignment is found, print it in alphabetical order
        result = list(perm)
        print(f"<<<{result}>>>")
        break