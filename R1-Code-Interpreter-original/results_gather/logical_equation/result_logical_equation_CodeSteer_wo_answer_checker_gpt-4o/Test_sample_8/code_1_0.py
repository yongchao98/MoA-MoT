from itertools import permutations

# The numbers to be assigned
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (B > F and
        D + E == 107 and
        D > F and
        D - F == 77 and
        C > B and
        C + I == 52 and
        B == 3 * F):
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I]
        print(f"<<<{result}>>>")
        break