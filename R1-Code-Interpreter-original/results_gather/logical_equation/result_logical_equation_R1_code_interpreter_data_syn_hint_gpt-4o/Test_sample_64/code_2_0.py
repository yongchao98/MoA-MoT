from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Iterate over possible values for L
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
    
    # Remaining numbers after assigning C, E, I, J, L
    remaining_numbers = [n for n in numbers if n not in {C, E, I, J, L}]
    
    # Check the remaining conditions
    for perm in permutations(remaining_numbers):
        A, B, D, F, G, H, K, M = perm
        
        if (E + H == 31 and
            F - K == 4 and
            D + K == 31 and
            D + G == 33 and
            A - H == -14 and
            J > G):
            
            # If all conditions are satisfied, print the result
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break