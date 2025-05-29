# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Function to check if a given assignment satisfies all conditions
def check_conditions(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        A - C == 26 and
        M > L and
        C == 2.0 * D and
        E + I == 57 and
        L == 1.6 * C and
        G - K == -12 and
        K == 1.5 * C and
        A == 3.6 * C and
        M == 4.8 * D and
        C + F == 106 and
        H - G == 42
    )

# Calculate possible values for C based on the equations
possible_C_values = [c for c in numbers if 3.6 * c in numbers and 2.0 * c / 2 in numbers and 1.5 * c in numbers and 1.6 * c in numbers]

# Try permutations of the remaining numbers for the other variables
from itertools import permutations

for C in possible_C_values:
    A = 3.6 * C
    D = C / 2
    K = 1.5 * C
    L = 1.6 * C
    M = 4.8 * D
    remaining_numbers = [n for n in numbers if n not in {A, C, D, K, L, M}]
    
    for perm in permutations(remaining_numbers):
        B, E, F, G, H, I, J = perm
        if check_conditions([A, B, C, D, E, F, G, H, I, J, K, L, M]):
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break