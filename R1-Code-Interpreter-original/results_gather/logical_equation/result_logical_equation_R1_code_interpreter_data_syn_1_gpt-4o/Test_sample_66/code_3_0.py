from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Pre-calculate possible values for D, I, M, L, G, A, and K
possible_DI = [(d, i) for d in numbers for i in numbers if d == 3 * i]
possible_ML = [(m, l) for m in numbers for l in numbers if m == 4.5 * l]
possible_LG = [(l, g) for l in numbers for g in numbers if l == g - 26]
possible_AK = [(a, k) for a in numbers for k in numbers if a == k + 17]

# Iterate over all possible combinations of pre-calculated values
found = False
for D, I in possible_DI:
    for M, L in possible_ML:
        for L2, G in possible_LG:
            if L != L2:
                continue
            for A, K in possible_AK:
                # Remaining numbers to permute
                remaining_numbers = set(numbers) - {A, D, G, I, K, L, M}
                
                # Check if the number of remaining numbers matches the expected number of variables
                if len(remaining_numbers) != 6:
                    continue
                
                # Iterate over permutations of the remaining numbers
                for perm in permutations(remaining_numbers):
                    B, C, E, F, H, J = perm
                    
                    # Check all conditions
                    if (C + H == 20 and
                        A > F and
                        F + G == 38 and
                        C + F == 17 and
                        E + L == 60 and
                        D + H == 8):
                        
                        # If all conditions are satisfied, print the result
                        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
                        found = True
                        break
                if found:
                    break
            if found:
                break
        if found:
            break
    if found:
        break