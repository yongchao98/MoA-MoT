from itertools import permutations

# The numbers to be assigned
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Iterate over all permutations of the numbers
for perm in permutations(numbers):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check if the current permutation satisfies all the inequalities
    if (K - G == 25 and
        F == 3.5 * M and
        B + E == 141 and
        J + K == 38 and
        B + I == 50 and
        A - G == 33 and
        I + M == 7 and
        A - H == -14 and
        M - H == -48 and
        C - E == -72 and
        F == 1.4 * I):
        
        # If all inequalities are satisfied, print the result
        result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
        print(f"<<<{result}>>>")
        break