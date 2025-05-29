from itertools import permutations

# Define the possible values for each letter
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Define the constraints as functions
def constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    return (
        J + K == 271 and
        G > I and
        E + J == 141 and
        A + C == 300 and
        G + K == 230 and
        B - K == -141 and
        D - C == -209 and
        A + I == 78 and
        F - D == 23 and
        E - H == -28 and
        H + K == 198
    )

# Try all permutations of the possible values
for perm in permutations(possible_values):
    if constraints(perm):
        solution = perm
        break

# Print the solution in alphabetical order
print(list(solution))