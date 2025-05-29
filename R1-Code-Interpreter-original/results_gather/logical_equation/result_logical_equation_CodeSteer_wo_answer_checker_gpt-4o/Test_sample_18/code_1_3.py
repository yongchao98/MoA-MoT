from itertools import permutations

# Define the numbers to be used
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check if the current permutation satisfies all the conditions
    if (A + C == 127 and
        D - E == -77 and
        A - B == 89 and
        B == 1.6 * F and
        E == 4.0 * F and
        I == 4.0 * H and
        H == 1.5 * C and
        E > I):
        
        # If all conditions are satisfied, print the result
        print(f"<<<{[A, B, C, D, E, F, G, H, I]}>>>")
        break