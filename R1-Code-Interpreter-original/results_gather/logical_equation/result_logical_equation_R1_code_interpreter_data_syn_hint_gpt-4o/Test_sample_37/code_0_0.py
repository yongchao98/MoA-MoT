from itertools import permutations

# List of numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if (I == 1.5 * H and
        F == 4.0 * B and
        C + E == 169 and
        E == 3.0 * J and
        J - G == -11 and
        C > F and
        I + K == 305 and
        I - D == 222 and
        J - F == -20 and
        F + H == 186):
        print([A, B, C, D, E, F, G, H, I, J, K])
        break