# Initial quantities
A = 6
B = 3
C = 4
X = 0

# Function to apply the methods
def apply_methods(A, B, C, X):
    # Method 1: 1A + 1B + 1C -> 1X
    if A >= 1 and B >= 1 and C >= 1:
        A -= 1
        B -= 1
        C -= 1
        X += 1
        return A, B, C, X, True
    
    # Method 2: 1A + 1B -> 2C
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        C += 2
        return A, B, C, X, True
    
    # Method 3: 2C -> 1X
    if C >= 2:
        C -= 2
        X += 1
        return A, B, C, X, True
    
    # No method can be applied
    return A, B, C, X, False

# Apply methods in cycle
while True:
    A, B, C, X, changed = apply_methods(A, B, C, X)
    if not changed:
        break

# Output the final result
print([A, B, C, X])