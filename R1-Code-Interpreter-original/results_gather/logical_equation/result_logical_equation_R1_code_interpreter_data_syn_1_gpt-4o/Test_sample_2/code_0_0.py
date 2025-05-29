from itertools import permutations

# The set of numbers to assign
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all the conditions
    if (G - E == -3 and
        D - C == -41 and
        C > B and
        F == 4 * H and
        I - G == 2 and
        A + E == 31 and
        A - B == 11):
        print([A, B, C, D, E, F, G, H, I])
        break