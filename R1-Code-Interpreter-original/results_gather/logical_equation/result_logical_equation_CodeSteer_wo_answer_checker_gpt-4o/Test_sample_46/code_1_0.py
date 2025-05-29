from itertools import permutations

# Define the numbers and the letters
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

# Define the constraints as a function
def satisfies_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    return (
        A + F == 196 and
        G - D == -6 and
        G - A == -72 and
        E == 3.0 * A and
        B > I and
        K - A == 75 and
        H - F == -41 and
        B - G == 36 and
        A - H == -5 and
        E - A == 150 and
        H - J == 64
    )

# Try all permutations of the numbers
for perm in permutations(numbers):
    if satisfies_constraints(perm):
        result = list(perm)
        break

# Output the result
print(f"<<<{result}>>>")