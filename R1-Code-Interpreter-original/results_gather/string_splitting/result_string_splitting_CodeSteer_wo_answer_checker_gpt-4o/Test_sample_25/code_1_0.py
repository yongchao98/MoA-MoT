# Initial quantities
A, B, C = 3, 1, 3
X, Y, Z = 0, 0, 0

# Function to check if any method can be executed
def can_execute_any_method(A, B, C, X, Y):
    return (A >= 1) or (B >= 2) or (C >= 2) or (B >= 1 and C >= 1) or (X >= 1 and Y >= 1)

# Loop until no method can be executed
while can_execute_any_method(A, B, C, X, Y):
    # Method 1: One unit of machine A -> 2 units of part X and 1 unit of part Y
    if A >= 1:
        A -= 1
        X += 2
        Y += 1
    # Method 2: Two units of machine B -> 1 unit of part X
    elif B >= 2:
        B -= 2
        X += 1
    # Method 3: Two units of machine C -> 1 unit of part Y
    elif C >= 2:
        C -= 2
        Y += 1
    # Method 4: One unit of machine B and one unit of machine C -> 1 unit of machine A
    elif B >= 1 and C >= 1:
        B -= 1
        C -= 1
        A += 1
    # Method 5: One unit of part X and one unit of part Y -> 1 unit of part Z
    elif X >= 1 and Y >= 1:
        X -= 1
        Y -= 1
        Z += 1

# Output the final quantities
result = [str(A), str(B), str(C), str(X), str(Y), str(Z)]
print(result)