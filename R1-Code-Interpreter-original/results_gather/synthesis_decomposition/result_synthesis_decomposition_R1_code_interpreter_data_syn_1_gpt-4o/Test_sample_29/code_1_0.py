# Initial quantities
A = 6
B = 4
X = 0

# Method execution loop
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
    else:
        continue  # Skip to next cycle if Method 2 cannot be executed

# Final result
print([A, B, X])