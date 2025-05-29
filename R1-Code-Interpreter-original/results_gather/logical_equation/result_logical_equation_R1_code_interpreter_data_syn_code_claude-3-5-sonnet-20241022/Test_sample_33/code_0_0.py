from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    # Assign values to letters
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (E - F == -59 and
        G - C == 6 and
        H + J == 107 and
        E + I == 137 and
        B == 4.0 * G and
        A + G == 234 and
        G == 3.0 * C and
        E - H == -64 and
        C - E == -13 and
        F + K == 123):
        
        # Print result in alphabetical order
        result = [A, B, C, D, E, F, G, H, I, J, K]
        print(result)
        break