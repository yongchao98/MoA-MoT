# Initial counts
A, B, C = 3, 3, 3
X, Y, Z = 0, 0, 0

# Function to simulate the process
def simulate_dismantling(A, B, C, X, Y, Z):
    while True:
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
        else:
            break
    return A, B, C, X, Y, Z

# Simulate the process
final_A, final_B, final_C, final_X, final_Y, final_Z = simulate_dismantling(A, B, C, X, Y, Z)

# Output the final result
print([str(final_A), str(final_B), str(final_C), str(final_X), str(final_Y), str(final_Z)])