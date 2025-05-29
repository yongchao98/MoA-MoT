def dismantle_machines():
    # Initial counts
    A, B, C = 4, 2, 4
    X, Y, Z = 0, 0, 0

    # Flag to check if any method was executed in the last cycle
    any_method_executed = True

    while any_method_executed:
        any_method_executed = False

        # Method 1
        if A >= 1:
            A -= 1
            X += 2
            Y += 1
            any_method_executed = True

        # Method 2
        if B >= 2:
            B -= 2
            X += 1
            any_method_executed = True

        # Method 3
        if C >= 2:
            C -= 2
            Y += 1
            any_method_executed = True

        # Method 4
        if B >= 1 and C >= 1:
            B -= 1
            C -= 1
            A += 1
            any_method_executed = True

        # Method 5
        if X >= 1 and Y >= 1:
            X -= 1
            Y -= 1
            Z += 1
            any_method_executed = True

    # Output the final counts
    result = [str(A), str(B), str(C), str(X), str(Y), str(Z)]
    print(result)

dismantle_machines()