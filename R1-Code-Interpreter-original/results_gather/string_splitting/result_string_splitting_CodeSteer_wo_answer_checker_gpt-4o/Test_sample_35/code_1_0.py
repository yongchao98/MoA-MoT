def dismantle_machines():
    # Initial counts
    A, B, C = 3, 3, 3
    X, Y, Z = 0, 0, 0

    # Loop until no method can be executed
    while True:
        executed = False

        # Method 1
        if A >= 1:
            A -= 1
            X += 2
            Y += 1
            executed = True

        # Method 2
        elif B >= 2:
            B -= 2
            X += 1
            executed = True

        # Method 3
        elif C >= 2:
            C -= 2
            Y += 1
            executed = True

        # Method 4
        elif B >= 1 and C >= 1:
            B -= 1
            C -= 1
            A += 1
            executed = True

        # Method 5
        elif X >= 1 and Y >= 1:
            X -= 1
            Y -= 1
            Z += 1
            executed = True

        # If no method was executed, break the loop
        if not executed:
            break

    # Return the final counts
    return [str(A), str(B), str(C), str(X), str(Y), str(Z)]

# Get the result and print it
result = dismantle_machines()
print(result)