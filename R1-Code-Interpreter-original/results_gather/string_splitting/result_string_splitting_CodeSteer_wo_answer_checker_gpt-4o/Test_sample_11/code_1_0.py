def dismantle_machines():
    # Initial quantities
    A, B, C = 3, 3, 1
    X, Y, Z = 0, 0, 0

    # Define the methods
    def method1():
        nonlocal A, X, Y
        if A >= 1:
            A -= 1
            X += 2
            Y += 1
            return True
        return False

    def method2():
        nonlocal B, X
        if B >= 2:
            B -= 2
            X += 1
            return True
        return False

    def method3():
        nonlocal C, Y
        if C >= 2:
            C -= 2
            Y += 1
            return True
        return False

    def method4():
        nonlocal A, B, C
        if B >= 1 and C >= 1:
            B -= 1
            C -= 1
            A += 1
            return True
        return False

    def method5():
        nonlocal X, Y, Z
        if X >= 1 and Y >= 1:
            X -= 1
            Y -= 1
            Z += 1
            return True
        return False

    # Cycle through the methods
    while True:
        executed = False
        for method in [method1, method2, method3, method4, method5]:
            if method():
                executed = True
        if not executed:
            break

    # Format the result
    result = [str(A), str(B), str(C), str(X), str(Y), str(Z)]
    print(result)

dismantle_machines()