from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Constraints as functions
def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    return (
        J - B == 27 and
        K - H == 123 and
        C - I == -71 and
        I - C == 71 and
        K == 2 * J and
        C - B == -39 and
        G - E == -33 and
        A + D == 346 and
        C == 3 * G and
        E - B == -12
    )

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        result = list(perm)
        break

# Print the result
print(f"<<<{result}>>>")