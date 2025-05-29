from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (D > A and D - A == 66 and F - G == -145 and B - G == -209 and
        E == 3.0 * B and I - A == 112 and D > C and A + H == 159 and
        D + H == 225 and E == 2.4 * J and A - H == -141):
        print([A, B, C, D, E, F, G, H, I, J, K])
        break