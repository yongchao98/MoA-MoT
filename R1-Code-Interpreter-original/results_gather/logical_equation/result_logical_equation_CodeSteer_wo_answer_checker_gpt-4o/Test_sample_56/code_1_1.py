from itertools import permutations

# Define the numbers and letters
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

# Define a function to check if the current assignment satisfies all constraints
def satisfies_constraints(assignment):
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

# Try all permutations of the numbers
for perm in permutations(numbers):
    if satisfies_constraints(perm):
        # If a valid assignment is found, print it
        print(f"<<<{list(perm)}>>>")
        break