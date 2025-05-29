from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    # Assign numbers to letters
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (H + J == 230 and
        C + I == 137 and
        C + G == 91 and
        G == 3.0 * F and
        B - E == -27 and
        J - C == 134 and
        E + G == 123 and
        D + F == 75 and
        E + K == 261 and
        H + K == 305):
        
        # Print result in required format
        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}]>>>")
        break