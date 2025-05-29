from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Possible values for A based on B = 4.5A and A - M = 7
possible_A = [x for x in numbers if 4.5 * x in numbers and x - 7 in numbers]

# Iterate over possible values for A
for A in possible_A:
    B = 4.5 * A
    M = A - 7
    if B not in numbers or M not in numbers:
        continue
    
    # Possible values for K based on A + K = 15
    K = 15 - A
    if K not in numbers:
        continue
    
    # Possible values for L based on M = 3.0L
    L = M / 3.0
    if L not in numbers:
        continue
    
    # Remaining numbers after assigning A, B, K, L, M
    remaining_numbers = [x for x in numbers if x not in {A, B, K, L, M}]
    
    # Iterate over permutations of the remaining numbers
    for perm in permutations(remaining_numbers):
        C, D, E, F, G, H, I, J = perm
        
        # Check remaining conditions
        if (F + J == 86 and
            H == 1.5 * A and
            J > D and
            B + E == 52 and
            A > G and
            H > L):
            # If all conditions are satisfied, print the result
            print([A, B, C, D, E, F, G, H, I, J, K, L, M])
            break