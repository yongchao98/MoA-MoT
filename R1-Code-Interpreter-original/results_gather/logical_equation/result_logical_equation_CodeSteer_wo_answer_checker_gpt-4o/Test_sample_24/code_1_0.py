from itertools import permutations

# Given numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (G == 1.5 * H and
        E - D == -41 and
        F == 1.5 * E and
        D + G == 130 and
        F + H == 126 and
        F + I == 123 and
        E == 2.5 * B and
        D - E == 41):
        
        # If all conditions are satisfied, print the result
        print(f"<<<{[A, B, C, D, E, F, G, H, I]}>>>")
        break