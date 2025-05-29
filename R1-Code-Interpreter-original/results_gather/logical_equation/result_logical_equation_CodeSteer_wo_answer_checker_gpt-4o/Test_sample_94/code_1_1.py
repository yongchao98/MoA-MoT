from itertools import permutations

# The numbers to be assigned
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if (G + I == 55 and
        E > D and
        G + L == 57 and
        G - E == 5 and
        C + H == 27 and
        A - I == 5 and
        F - I == 91 and
        B - C == -9 and
        K > C and
        H + L == 10 and
        J - D == 14 and
        B - F == -81):
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
        print(f"<<<{result}>>>")
        break