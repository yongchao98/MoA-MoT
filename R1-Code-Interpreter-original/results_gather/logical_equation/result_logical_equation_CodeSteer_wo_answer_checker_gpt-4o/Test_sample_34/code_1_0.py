from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (G == 3.0 * C and
        I == 4.0 * G and
        A - G == 112 and
        K == 3.0 * J and
        B + D == 64 and
        E - B == 11 and
        F == 2.0 * J and
        F - C == 147 and
        E - I == -9 and
        A + F == 271):
        
        # If all conditions are satisfied, print the result
        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}]>>>")
        break