# Initial counts
A = 6
B = 3
C = 3
X = 0

# Function to simulate the exchange process
def simulate_exchange(A, B, C, X):
    while True:
        # Method 1: 1A + 1B + 1C -> 1X
        if A >= 1 and B >= 1 and C >= 1:
            A -= 1
            B -= 1
            C -= 1
            X += 1
        # Method 2: 1A + 1B -> 2C
        elif A >= 1 and B >= 1:
            A -= 1
            B -= 1
            C += 2
        # Method 3: 2C -> 1X
        elif C >= 2:
            C -= 2
            X += 1
        else:
            # No more exchanges can be made
            break
    return A, B, C, X

# Run the simulation
remaining_A, remaining_B, remaining_C, obtained_X = simulate_exchange(A, B, C, X)

# Output the result
print([remaining_A, remaining_B, remaining_C, obtained_X])