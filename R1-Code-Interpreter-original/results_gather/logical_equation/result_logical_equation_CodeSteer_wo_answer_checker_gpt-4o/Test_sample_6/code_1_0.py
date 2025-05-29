from itertools import permutations

# Possible values for the letters
values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (E + F == 83 and
        A - H == -112 and
        G == 4 * C and
        B + D == 63 and
        D == 3 * A and
        G + H == 125 and
        B - E == 33):
        
        # If all conditions are satisfied, print the result
        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}]>>>")
        break