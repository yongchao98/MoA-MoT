from itertools import permutations

# Define the possible values
values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Iterate over all permutations of the values
for perm in permutations(values):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check if the current permutation satisfies all the equations
    if (K - H == -150 and
        G - K == 46 and
        B - F == 70 and
        K - A == 59 and
        H == 3 * K and
        C + H == 252 and
        D == 4 * I and
        C + F == 107 and
        J - B == -102 and
        G - I == 112):
        
        # If all conditions are satisfied, print the solution
        solution = [A, B, C, D, E, F, G, H, I, J, K]
        break

print(solution)