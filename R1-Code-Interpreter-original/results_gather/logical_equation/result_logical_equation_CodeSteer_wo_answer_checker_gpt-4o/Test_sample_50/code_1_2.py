from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all the conditions
    if (B + F == 346 and
        A + E == 155 and
        F + G == 228 and
        A + H == 128 and
        A == 4 * D and
        I - K == -23 and
        J == 2 * E and
        F == 3 * E and
        B - J == -29 and
        A + K == 119 and
        H == 3 * I):
        
        # If all conditions are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K]
        print(f"<<<{result}>>>")
        break