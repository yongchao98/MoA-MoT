from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    # Assign values to letters
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (C - K == -17 and
        J == 2.0 * D and
        G - K == 205 and
        E + H == 119 and
        C - J == -147 and
        F + K == 36 and
        F + G == 241 and
        H - F == 64 and
        I - B == -39 and
        D + K == 95 and
        B + J == 198):
        
        # Print result in required format
        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}]>>>")
        break