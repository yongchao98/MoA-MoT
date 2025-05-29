from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    # Assign numbers to letters
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (H + J == 130 and
        C + I == 137 and
        C + G == 75 and
        abs(G - 3.0 * F) < 0.0001 and  # Using approximate equality for floating point
        B - E == 150 and
        J - C == -18 and
        E + G == 123 and
        D + F == 96 and
        E + K == 78 and
        H + K == 124):
        
        # Print result in required format
        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}]>>>")
        break