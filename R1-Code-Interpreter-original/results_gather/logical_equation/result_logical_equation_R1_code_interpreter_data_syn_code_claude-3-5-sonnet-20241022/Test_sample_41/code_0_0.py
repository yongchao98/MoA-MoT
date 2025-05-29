from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    # Assign values to letters
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all equations
    if (
        abs(C - 3.0*J) < 0.1 and           # C = 3.0J
        F + K == 375 and                    # F + K = 375
        abs(H - 3.0*D) < 0.1 and           # H = 3.0D
        D + F == 234 and                    # D + F = 234
        J - B == -59 and                    # J - B = -59
        J + K == 166 and                    # J + K = 166
        B - H == 48 and                     # B - H = 48
        K - C == 102 and                    # K - C = 102
        C + K == 198 and                    # C + K = 198
        D - B == -66                        # D - B = -66
    ):
        # Print result in alphabetical order
        print([A, B, C, D, E, F, G, H, I, J, K])
        break