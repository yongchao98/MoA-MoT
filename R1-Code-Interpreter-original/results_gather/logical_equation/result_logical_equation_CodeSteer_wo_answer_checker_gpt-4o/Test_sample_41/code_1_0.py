from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Define the equations as functions
def check_equations(values):
    A, B, C, D, E, F, G, H, I, J, K = values
    return (
        C == 3 * J and
        F + K == 375 and
        H == 3 * D and
        D + F == 234 and
        J - B == -59 and
        J + K == 166 and
        B - H == 48 and
        K - C == 102 and
        C + K == 198 and
        D - B == -66
    )

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    if check_equations(perm):
        # Output the results in alphabetical order
        print(f"<<<{list(perm)}>>>")
        break