# Initial quantities
A = 5
B = 3
C = 5
X = 0

# Function to apply the methods in order
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
        return A, B, C, X, False  # No more methods can be applied
    return A, B, C, X, True

# Apply methods until no more can be applied
while True:
    A, B, C, X, can_continue = apply_methods(A, B, C, X)
    if not can_continue:
        break

# Output the final state
print([A, B, C, X])