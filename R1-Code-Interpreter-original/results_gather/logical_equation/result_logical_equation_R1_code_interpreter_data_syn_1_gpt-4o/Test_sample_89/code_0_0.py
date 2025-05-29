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

# Try all permutations of the numbers to find a valid assignment
from itertools import permutations

for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break