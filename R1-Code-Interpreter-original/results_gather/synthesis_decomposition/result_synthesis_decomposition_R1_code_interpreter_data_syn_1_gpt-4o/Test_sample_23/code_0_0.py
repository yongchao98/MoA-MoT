# Initial quantities
A = 4
B = 4
X = 0

# Method application cycle
while True:
    # Try Method 1
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        X += 1
    else:
        break  # Cannot apply Method 1, exit the loop

    # Try Method 2
    if A >= 3:
        A -= 3
        X += 2
    else:
        break  # Cannot apply Method 2, exit the loop

# Output the final state
print([A, B, X])