# Initial quantities
A = 8
B = 6
C = 4
X = 0

# Function to perform the cycle of exchanges
def perform_exchanges(A, B, C, X):
    while True:
        # Method 1
        if A >= 1 and B >= 1 and C >= 1:
            A -= 1
            B -= 1
            C -= 1
            X += 1
        else:
            break
        
        # Method 2
        if A >= 1 and B >= 1:
            A -= 1
            B -= 1
            C += 2
        else:
            break
        
        # Method 3
        if C >= 2:
            C -= 2
            X += 1
        else:
            break
    
    return A, B, C, X

# Perform the exchanges
remaining_A, remaining_B, remaining_C, obtained_X = perform_exchanges(A, B, C, X)

# Output the result
print([remaining_A, remaining_B, remaining_C, obtained_X])