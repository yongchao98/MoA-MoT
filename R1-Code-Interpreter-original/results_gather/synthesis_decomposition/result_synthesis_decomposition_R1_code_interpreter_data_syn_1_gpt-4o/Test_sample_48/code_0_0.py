# Initial quantities
A = 8
B = 5
C = 4
X = 0

# Function to apply the methods
def apply_methods(A, B, C, X):
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
        return A, B, C, X, False  # No method can be applied
    return A, B, C, X, True  # A method was applied

# Simulate the process
while True:
    A, B, C, X, method_applied = apply_methods(A, B, C, X)
    if not method_applied:
        break

# Output the final state
print([A, B, C, X])