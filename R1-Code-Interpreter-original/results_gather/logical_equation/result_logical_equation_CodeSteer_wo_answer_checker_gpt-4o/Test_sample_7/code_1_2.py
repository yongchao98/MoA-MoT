from itertools import permutations

# Define the set of possible numbers
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check if the current permutation satisfies all the constraints
    if (D + H == 13 and
        E - C == 15 and
        C + H == 5 and
        H == 4 * C and
        H - B == -117 and
        F == 3 * C and
        F + H == 7):
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I]
        print(f"<<<{result}>>>")
        break