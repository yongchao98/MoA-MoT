from itertools import permutations

# The numbers to assign
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if (B == 2.0 * G and
        I == 2.5 * F and
        B + I == 86 and
        A + G == 12 and
        I > F and
        F + I == 112 and
        C + H == 140 and
        A - H == -11):
        
        # If all constraints are satisfied, print the solution
        print([A, B, C, D, E, F, G, H, I])
        break