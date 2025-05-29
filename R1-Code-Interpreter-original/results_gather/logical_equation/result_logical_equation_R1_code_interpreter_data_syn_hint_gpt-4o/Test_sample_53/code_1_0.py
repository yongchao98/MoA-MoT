from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all the conditions
    if (F - K == 36 and
        B - K == 111 and
        B + E == 230 and
        H + J == 25 and
        E - A == 77 and
        H + I == 36 and
        C + E == 305 and
        C - I == 205 and
        C + J == 234 and
        G == 3 * H and
        A + K == 42):
        
        # If all conditions are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K]
        print(result)
        break