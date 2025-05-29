from itertools import permutations

# The numbers to be assigned
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check the constraints
    if (C - D == -94 and
        C - F == 18 and
        B > E and
        A > E and
        G == 4 * E and
        B > F and
        C == 3 * F):
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I]
        print(f"<<<{result}>>>")
        break