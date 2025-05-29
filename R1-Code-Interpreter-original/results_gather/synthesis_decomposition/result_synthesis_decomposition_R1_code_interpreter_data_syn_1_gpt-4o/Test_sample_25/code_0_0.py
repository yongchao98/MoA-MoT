# Initial quantities
A = 3
B = 3
X = 0

# Method application loop
while True:
    # Try Method 1
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        X += 1
    else:
        break  # Can't apply Method 1, break the loop

    # Try Method 2
    if A >= 3:
        A -= 3
        X += 2
    else:
        continue  # Can't apply Method 2, continue to next cycle

# Output the final state
print([A, B, X])