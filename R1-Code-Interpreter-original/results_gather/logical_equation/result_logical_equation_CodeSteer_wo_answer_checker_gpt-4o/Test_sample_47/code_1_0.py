from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Function to check if a given permutation satisfies all constraints
def satisfies_constraints(perm):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    return (
        I == 4 * C and
        F + I == 83 and
        K - E == 82 and
        E + F == 42 and
        A + H == 225 and
        B + K == 130 and
        A > B and
        J == 3 * H and
        B == 3 * F and
        C + F == 23 and
        A > K
    )

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    if satisfies_constraints(perm):
        # If a valid permutation is found, print the result
        print(f"<<<{list(perm)}>>>")
        break