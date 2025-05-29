def optimize_dismantling():
    # Initial quantities
    A, B, C = 4, 4, 4
    X, Y, Z = 0, 0, 0

    # Function to evaluate and execute the best method
    def execute_best_method():
        nonlocal A, B, C, X, Y, Z
        # List of possible actions with their conditions and effects
        actions = [
            ("Method 1", A >= 1, lambda: (A - 1, B, C, X + 2, Y + 1, Z)),
            ("Method 2", B >= 2, lambda: (A, B - 2, C, X + 1, Y, Z)),
            ("Method 3", C >= 2, lambda: (A, B, C - 2, X, Y + 1, Z)),
            ("Method 4", B >= 1 and C >= 1, lambda: (A + 1, B - 1, C - 1, X, Y, Z)),
            ("Method 5", X >= 1 and Y >= 1, lambda: (A, B, C, X - 1, Y - 1, Z + 1))
        ]

        # Evaluate and execute the best action
        for name, condition, effect in actions:
            if condition:
                A, B, C, X, Y, Z = effect()
                return True
        return False

    # Loop until no method can be executed
    while execute_best_method():
        pass

    # Output the final quantities
    result = [str(A), str(B), str(C), str(X), str(Y), str(Z)]
    print(result)

optimize_dismantling()