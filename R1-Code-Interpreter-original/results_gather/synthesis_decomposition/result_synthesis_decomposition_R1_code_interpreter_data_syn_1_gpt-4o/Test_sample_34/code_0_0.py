# Initial quantities
A = 8
B = 5
C = 2
X = 0

# Function to execute the methods
def execute_methods(A, B, C, X):
    while True:
        # Method 1
        if A >= 1 and B >= 1 and C >= 1:
            A -= 1
            B -= 1
            C -= 1
            X += 1
        # Method 2
        elif A >= 1 and B >= 1:
            A -= 1
            B -= 1
            C += 2
        # Method 3
        elif C >= 2:
            C -= 2
            X += 1
        else:
            break
    return A, B, C, X

# Execute the methods
remaining_A, remaining_B, remaining_C, obtained_X = execute_methods(A, B, C, X)

# Output the final result
print([remaining_A, remaining_B, remaining_C, obtained_X])