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
        return A, B, C, X, True
    
    # Method 2
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        C += 2
        return A, B, C, X, True
    
    # Method 3
    if C >= 2:
        C -= 2
        X += 1
        return A, B, C, X, True
    
    # If no method can be applied
    return A, B, C, X, False

# Apply methods in a cycle
while True:
    A, B, C, X, changed = apply_methods(A, B, C, X)
    if not changed:
        break

# Output the final state
print([A, B, C, X])