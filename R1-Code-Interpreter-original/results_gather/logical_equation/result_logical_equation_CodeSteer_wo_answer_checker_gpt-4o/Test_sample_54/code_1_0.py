import itertools

# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in itertools.permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all constraints
    if (H > G and
        H + K == 273 and
        B - F == 60 and
        K - B == 145 and
        F - I == 11 and
        D - E == -59 and
        B == 4 * F and
        G + H == 51 and
        K == 3 * E and
        C + K == 264 and
        I == 3 * G):
        
        # If all constraints are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K]
        print(f"<<<{result}>>>")
        break