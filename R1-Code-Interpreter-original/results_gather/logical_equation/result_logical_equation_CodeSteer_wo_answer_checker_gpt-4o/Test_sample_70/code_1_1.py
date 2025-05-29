from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if (M == 2.4 * L and
        F == 3.0 * L and
        E + I == 8 and
        K - H == -8 and
        H > I and
        L + M == 51 and
        C + J == 9 and
        G == 2.0 * I and
        A - D == 49 and
        H == 1.5 * K):
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
        print(f"<<<{result}>>>")
        break