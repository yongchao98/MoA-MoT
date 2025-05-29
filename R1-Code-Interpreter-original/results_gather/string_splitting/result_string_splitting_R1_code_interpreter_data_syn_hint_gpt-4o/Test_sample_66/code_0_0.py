# Initial quantities
A, B, C = 3, 2, 2
X, Y, Z = 0, 0, 0

# Function to check if any method can be executed
def can_execute_any_method(A, B, C, X, Y):
    return (A > 0) or (B >= 2) or (C >= 2) or (B > 0 and C > 0) or (X > 0 and Y > 0)

# Simulate the process
while can_execute_any_method(A, B, C, X, Y):
    # Method 1
    if A > 0:
        A -= 1
        X += 2
        Y += 1
    # Method 2
    elif B >= 2:
        B -= 2
        X += 1
    # Method 3
    elif C >= 2:
        C -= 2
        Y += 1
    # Method 4
    elif B > 0 and C > 0:
        B -= 1
        C -= 1
        A += 1
    # Method 5
    elif X > 0 and Y > 0:
        X -= 1
        Y -= 1
        Z += 1

# Final quantities
result = [str(A), str(B), str(C), str(X), str(Y), str(Z)]
print(result)