# Initial quantities
A = 5
B = 6
C = 4
X = 0

# Function to execute the methods in cycle
def execute_methods(A, B, C, X):
    while True:
        # Method 1: 1A + 1B + 1C -> 1X
        if A >= 1 and B >= 1 and C >= 1:
            A -= 1
            B -= 1
            C -= 1
            X += 1
        # Method 2: 1A + 1B -> 2C
        elif A >= 1 and B >= 1:
            A -= 1
            B -= 1
            C += 2
        # Method 3: 2C -> 1X
        elif C >= 2:
            C -= 2
            X += 1
        else:
            # No more methods can be executed
            break
    return A, B, C, X

# Execute the methods
remaining_A, remaining_B, remaining_C, obtained_X = execute_methods(A, B, C, X)

# Output the final result
print([remaining_A, remaining_B, remaining_C, obtained_X])