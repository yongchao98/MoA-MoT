# Initial quantities
A = 6
B = 4
X = 0

# Cycle through methods
while True:
    # Method 1
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        X += 1
    else:
        break

    # Method 2
    if A >= 3:
        A -= 3
        X += 2
    else:
        break

# Output the final state
print([A, B, X])