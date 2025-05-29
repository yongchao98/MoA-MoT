from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    # Assign values to letters
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all constraints
    if (K - H == -150 and
        G - K == 46 and
        B - F == 70 and
        K - A == 59 and
        H == 3.0 * K and
        C + H == 252 and
        D == 4.0 * I and
        C + F == 107 and
        J - B == -102 and
        G - I == 112):
        
        # Print result in alphabetical order
        result = [A, B, C, D, E, F, G, H, I, J, K]
        print(result)
        break