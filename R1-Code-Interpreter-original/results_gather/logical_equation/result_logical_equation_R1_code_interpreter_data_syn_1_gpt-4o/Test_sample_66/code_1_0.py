from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Pre-calculate possible values for D, I, M, L, G, A, and K
possible_DI = [(d, i) for d in numbers for i in numbers if d == 3 * i]
possible_ML = [(m, l) for m in numbers for l in numbers if m == 4.5 * l]
possible_LG = [(l, g) for l in numbers for g in numbers if l == g - 26]
possible_AK = [(a, k) for a in numbers for k in numbers if a == k + 17]

# Iterate over all permutations of the remaining numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check if the current permutation satisfies the pre-calculated constraints
    if ((D, I) in possible_DI and
        (M, L) in possible_ML and
        (L, G) in possible_LG and
        (A, K) in possible_AK and
        C + H == 20 and
        A > F and
        F + G == 38 and
        C + F == 17 and
        E + L == 60 and
        D + H == 8):
        
        # If all conditions are satisfied, print the result
        print([A, B, C, D, E, F, G, H, I, J, K, L, M])
        break