from itertools import permutations

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if (D + E == 15 and
        A + F == 201 and
        D == 1.5 * E and
        C == 4.0 * D and
        H > G and
        E - I == -26 and
        E == 2.0 * G and
        C + I == 68):
        
        # If all conditions are satisfied, print the result
        print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}]>>>")
        break