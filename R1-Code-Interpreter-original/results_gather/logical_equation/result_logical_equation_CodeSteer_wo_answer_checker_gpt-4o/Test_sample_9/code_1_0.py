import itertools

# The numbers to be assigned
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the numbers
for perm in itertools.permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (D - E == 24 and
        A - F == 35 and
        D == 3 * B and
        I == 4 * F and
        B == 3 * E and
        E + H == 124 and
        C > E):
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I]
        print(f"<<<{result}>>>")
        break