from itertools import permutations

# Given numbers
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Function to check if a permutation satisfies all constraints
def satisfies_constraints(perm):
    A, B, C, D, E, F, G, H, I = perm
    return (
        G == 4 * A and
        C == 4 * D and
        B + E == 148 and
        A - E == -117 and
        D - F == -71 and
        C + H == 39
    )

# Iterate through all permutations
for perm in permutations(numbers):
    if satisfies_constraints(perm):
        # Print the result in the required format
        print(f"<<<{list(perm)}>>>")
        break