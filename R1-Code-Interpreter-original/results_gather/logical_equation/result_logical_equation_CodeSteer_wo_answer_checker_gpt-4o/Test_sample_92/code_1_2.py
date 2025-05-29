from itertools import permutations

# Define the possible values for each letter
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Define the constraints as functions
def check_constraints(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    return (
        A == 3.0 * L and
        A + M == 39 and
        L + M == 29 and
        B == 3.2 * L and
        A == 1.5 * I and
        C - B == 80 and
        E + K == 43 and
        C > D and
        E == 2.4 * A and
        C + D == 124 and
        H - A == -12 and
        F + G == 47
    )

# Iterate over all permutations of the possible values
for perm in permutations(possible_values):
    if check_constraints(perm):
        result = list(perm)
        print(f"<<<{result}>>>")
        break