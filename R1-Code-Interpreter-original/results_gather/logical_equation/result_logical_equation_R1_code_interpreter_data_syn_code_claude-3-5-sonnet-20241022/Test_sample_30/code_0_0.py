from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    # Assign numbers to letters
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (C + J == 186 and
        abs(H - 3.0 * F) < 0.1 and  # Using small threshold for float comparison
        C - B == -44 and
        E + H == 234 and
        E + I == 130 and
        I - J == -29 and
        D + E == 57 and
        K - B == -64 and
        C + D == 84 and
        B - G == 53):
        
        # Print result in required format
        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}]>>>")
        break