from itertools import permutations

# List of numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Function to check if a permutation satisfies all the conditions
def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K = perm
    return (
        J - B == 27 and
        K - H == 123 and
        C - I == -71 and
        I - C == 71 and
        K == 2.0 * J and
        C - B == -39 and
        G - E == -33 and
        A + D == 346 and
        C == 3.0 * G and
        E - B == -12
    )

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break