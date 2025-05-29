from itertools import permutations

# List of numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (D + J == 252 and
        D == 3 * K and
        H - A == 44 and
        J > E and
        C == 2 * K and
        D > A and
        A + F == 52 and
        D + E == 228 and
        B == 3 * F and
        E + I == 124):
        print([A, B, C, D, E, F, G, H, I, J, K])
        break