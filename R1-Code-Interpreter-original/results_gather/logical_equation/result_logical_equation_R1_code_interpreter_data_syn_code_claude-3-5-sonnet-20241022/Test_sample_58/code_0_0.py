from itertools import permutations

# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Function to check if a permutation satisfies all conditions
def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    conditions = [
        B + E == 170,
        H - C == -46,
        K == 3.0 * G,
        B + I == 230,
        C - H == 46,
        D == 3.0 * H,
        C - D == -104,
        H - J == 36,
        A == 2.4 * E,
        I - E == 60,
        C - K == 112
    ]
    
    return all(conditions)

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break