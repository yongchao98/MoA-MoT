# Initial quantities
A = 6
B = 3
X = 0

# Function to execute the cycle
def execute_cycle(A, B, X):
    # Method 1
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        X += 1
    
    # Method 2
    if A >= 3:
        A -= 3
        X += 2
    
    return A, B, X

# Execute cycles until no more methods can be executed
while True:
    new_A, new_B, new_X = execute_cycle(A, B, X)
    if new_A == A and new_B == B and new_X == X:
        break
    A, B, X = new_A, new_B, new_X

print([A, B, X])