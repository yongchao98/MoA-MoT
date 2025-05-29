from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over possible values for L and E
for L in numbers:
    E = 1.5 * L
    if E not in numbers:
        continue
    
    I = 2.4 * E
    if I not in numbers:
        continue
    
    J = 1.6 * E
    if J not in numbers:
        continue
    
    C = I - 35
    if C not in numbers:
        continue
    
    # Now check the remaining conditions
    for perm in permutations(numbers):
        A, B, D, F, G, H, K, M = perm
        
        if (E + H == 31 and
            F - K == 4 and
            D + K == 31 and
            D + G == 33 and
            A - H == -14 and
            J > G and
            A != E and B != E and D != E and F != E and G != E and H != E and K != E and M != E and
            A != I and B != I and D != I and F != I and G != I and H != I and K != I and M != I and
            A != J and B != J and D != J and F != J and G != J and H != J and K != J and M != J and
            A != C and B != C and D != C and F != C and G != C and H != C and K != C and M != C):
            
            # If all conditions are satisfied, print the result
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break