import itertools

# The numbers to be assigned
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the numbers
for perm in itertools.permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (E + G == 127 and
        H == 4 * D and
        A + D == 140 and
        B == 1.6 * D and
        A - G == -1 and
        G > B and
        C + D == 29 and
        G - C == 112):
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I]
        print(f"<<<{result}>>>")
        break