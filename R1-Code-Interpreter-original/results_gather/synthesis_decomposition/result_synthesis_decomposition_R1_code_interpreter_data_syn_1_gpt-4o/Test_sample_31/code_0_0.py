# Initial quantities
A = 8
B = 4
C = 2
X = 0

# Function to perform the exchange cycle
def exchange_cycle(A, B, C, X):
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

# Perform the cycle until no method can be executed
while True:
    new_A, new_B, new_C, new_X = exchange_cycle(A, B, C, X)
    # Check if any changes occurred
    if (new_A, new_B, new_C, new_X) == (A, B, C, X):
        break
    A, B, C, X = new_A, new_B, new_C, new_X

# Output the final state
print([A, B, C, X])