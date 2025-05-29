# Initial quantities
A = 7
B = 4
C = 5
X = 0

# Function to apply the methods
def apply_methods(A, B, C, X):
    while True:
        # Method 1: 1A + 1B + 1C -> 1X
        if A >= 1 and B >= 1 and C >= 1:
            A -= 1
            B -= 1
            C -= 1
            X += 1
        else:
            # Method 2: 1A + 1B -> 2C
            if A >= 1 and B >= 1:
                A -= 1
                B -= 1
                C += 2
            else:
                # Method 3: 2C -> 1X
                if C >= 2:
                    C -= 2
                    X += 1
                else:
                    # No more methods can be applied
                    break
    return A, B, C, X

# Apply the methods
remaining_A, remaining_B, remaining_C, obtained_X = apply_methods(A, B, C, X)

# Output the result
print([remaining_A, remaining_B, remaining_C, obtained_X])