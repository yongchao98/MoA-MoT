from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Solve some constraints analytically
# 1. D = 3.0I
possible_DI = [(d, i) for d in numbers for i in numbers if d == 3 * i]

# 2. M = 4.5L
possible_ML = [(m, l) for m in numbers for l in numbers if m == 4.5 * l]

# 3. L - G = -26
possible_LG = [(l, g) for l in numbers for g in numbers if l - g == -26]

# 4. F + G = 38
possible_FG = [(f, g) for f in numbers for g in numbers if f + g == 38]

# 5. C + F = 17
possible_CF = [(c, f) for c in numbers for f in numbers if c + f == 17]

# 6. E + L = 60
possible_EL = [(e, l) for e in numbers for l in numbers if e + l == 60]

# 7. A - K = -17
possible_AK = [(a, k) for a in numbers for k in numbers if a - k == -17]

# 8. C + H = 20
possible_CH = [(c, h) for c in numbers for h in numbers if c + h == 20]

# 9. D + H = 8
possible_DH = [(d, h) for d in numbers for h in numbers if d + h == 8]

# Function to check if a partial assignment satisfies the constraints
def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        (D, I) in possible_DI and
        (M, L) in possible_ML and
        (L, G) in possible_LG and
        (F, G) in possible_FG and
        (C, F) in possible_CF and
        (E, L) in possible_EL and
        (A, K) in possible_AK and
        (C, H) in possible_CH and
        (D, H) in possible_DH and
        A > F
    )

# Backtracking function with constraint propagation
def backtrack(assignment, index):
    if index == len(numbers):
        if check_constraints(assignment):
            return assignment
        return None
    
    for num in numbers:
        if num not in assignment:
            assignment[index] = num
            result = backtrack(assignment, index + 1)
            if result is not None:
                return result
            assignment[index] = None
    
    return None

# Initialize assignment
assignment = [None] * 13

# Find solution using backtracking
solution = backtrack(assignment, 0)

# Output the solution
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")