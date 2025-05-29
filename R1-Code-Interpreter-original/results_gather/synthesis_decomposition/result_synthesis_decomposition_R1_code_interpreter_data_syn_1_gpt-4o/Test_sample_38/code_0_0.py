# Initial quantities
A = 8
B = 5
C = 5
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
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        C += 2
    
    # Method 3
    if C >= 2:
        C -= 2
        X += 1
    
    return A, B, C, X

# Apply the methods in cycles
while True:
    initial_state = (A, B, C, X)
    A, B, C, X = apply_methods(A, B, C, X)
    # If no changes occur, break the loop
    if initial_state == (A, B, C, X):
        break

# Output the final state
print([A, B, C, X])