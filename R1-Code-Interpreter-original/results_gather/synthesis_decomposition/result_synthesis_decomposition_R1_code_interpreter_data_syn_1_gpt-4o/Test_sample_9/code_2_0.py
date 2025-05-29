# Initial quantities
A = 3
B = 2
X = 0

# Method 1: 1A + 1B -> 1X
# Method 2: 3A -> 2X

# Execute the cycle
while True:
    # Try Method 1
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        X += 1
    else:
        break  # Cannot execute Method 1

    # Try Method 2
    if A >= 3:
        A -= 3
        X += 2

# Output the final state
print([A, B, X])