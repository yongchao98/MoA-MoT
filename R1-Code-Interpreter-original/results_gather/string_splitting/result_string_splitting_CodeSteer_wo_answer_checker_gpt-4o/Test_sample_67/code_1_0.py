# Initial quantities
A, B, C = 3, 2, 3
X, Y, Z = 0, 0, 0

# Function to check if any method can be executed
def can_execute():
    return (A >= 1 or B >= 2 or C >= 2 or (B >= 1 and C >= 1) or (X >= 1 and Y >= 1))

# Loop through the methods until no method can be executed
while can_execute():
    # Method 1
    if A >= 1:
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
    elif B >= 1 and C >= 1:
        B -= 1
        C -= 1
        A += 1
    # Method 5
    elif X >= 1 and Y >= 1:
        X -= 1
        Y -= 1
        Z += 1

# Output the final quantities
result = [str(A), str(B), str(C), str(X), str(Y), str(Z)]
print(result)