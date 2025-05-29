# Initial quantities
A = 7
B = 3
C = 2
X = 0

# Function to execute the methods
def execute_methods(A, B, C, X):
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
    
    # If no method can be executed
    return A, B, C, X, False

# Simulate the process
while True:
    A, B, C, X, executed = execute_methods(A, B, C, X)
    if not executed:
        break

# Output the final result
print([A, B, C, X])