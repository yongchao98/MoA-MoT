from itertools import permutations

# The numbers to be assigned
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, F, G = perm[:4]  # Choose 4 numbers for A, B, F, G
    
    # Calculate C and E based on the constraint C - E = -11
    for C in numbers:
        E = C + 11
        if E not in numbers:
            continue
        
        # Calculate J based on J = 4.8C
        J = 4.8 * C
        if J not in numbers:
            continue
        
        # Calculate D based on J - D = 14
        D = J - 14
        if D not in numbers:
            continue
        
        # Calculate H based on J - H = 22
        H = J - 22
        if H not in numbers:
            continue
        
        # Calculate L and then I and M based on I = L - 12 and M = 2.4L
        for L in numbers:
            I = L - 12
            M = 2.4 * L
            if I not in numbers or M not in numbers:
                continue
            
            # Calculate K based on K = 3.5H
            K = 3.5 * H
            if K not in numbers:
                continue
            
            # Check if all values are distinct
            if len(set([A, B, C, D, E, F, G, H, I, J, K, L, M])) == 13:
                # If all conditions are satisfied, print the result
                print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                break