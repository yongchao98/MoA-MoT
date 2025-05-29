from itertools import permutations

# Given numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (C == 3 * H and
        A - I == 16 and
        I - H == 17 and
        B == 4 * I and
        E + F == 127 and
        F - G == -114 and
        A == 4 * C and
        B == 2.5 * D):
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I]
        print(f"<<<{result}>>>")
        break