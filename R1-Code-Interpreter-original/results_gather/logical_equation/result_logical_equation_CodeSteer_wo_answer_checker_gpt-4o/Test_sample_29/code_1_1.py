from itertools import permutations

# Available numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, F, G, H, I = perm
    
    # Check all conditions
    if (B + I == 156 and
        B + D == 123 and
        D - A == -29 and
        G == 2 * D and
        B == 1.5 * C and
        A > D and
        A - I == -4 and
        C - F == 71):
        
        # Print the result in alphabetical order
        result = [A, B, C, D, F, G, H, I]
        print(f"<<<{result}>>>")
        break