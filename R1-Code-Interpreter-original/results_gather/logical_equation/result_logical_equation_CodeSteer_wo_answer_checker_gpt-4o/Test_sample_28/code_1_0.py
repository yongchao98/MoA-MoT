from itertools import permutations

possible_values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Iterate over all permutations of the possible values
for perm in permutations(possible_values):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check the constraints
    if (G == 1.5 * C and
        E == 1.8 * D and
        D + I == 141 and
        B - H == -88 and
        H + I == 241 and
        H - G == 111 and
        H - B == 88 and
        A + I == 201):
        
        print(f"<<<{[A, B, C, D, E, F, G, H, I]}>>>")
        break