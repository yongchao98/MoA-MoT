def dismantle_machines():
    # Initial counts
    A, B, C = 1, 3, 3
    X, Y, Z = 0, 0, 0

    while True:
        # Track if any method was executed in this cycle
        executed = False

        # Method 1: A -> 2X + Y
        if A >= 1:
            A -= 1
            X += 2
            Y += 1
            executed = True

        # Method 2: 2B -> X
        elif B >= 2:
            B -= 2
            X += 1
            executed = True

        # Method 3: 2C -> Y
        elif C >= 2:
            C -= 2
            Y += 1
            executed = True

        # Method 4: B + C -> A
        elif B >= 1 and C >= 1:
            B -= 1
            C -= 1
            A += 1
            executed = True

        # Method 5: X + Y -> Z
        elif X >= 1 and Y >= 1:
            X -= 1
            Y -= 1
            Z += 1
            executed = True

        # If no method was executed, break the loop
        if not executed:
            break

    # Return the final counts in the specified format
    return [str(A), str(B), str(C), str(X), str(Y), str(Z)]

# Execute the function and print the result
result = dismantle_machines()
print(result)