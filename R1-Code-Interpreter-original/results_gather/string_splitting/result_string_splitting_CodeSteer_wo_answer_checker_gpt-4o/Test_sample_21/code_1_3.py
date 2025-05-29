from sympy import symbols, Eq, solve, Min

def symbolic_dismantle_and_combine():
    # Define symbolic variables
    A, B, C, X, Y, Z = symbols('A B C X Y Z', integer=True)

    # Initial counts
    initial_A = 3
    initial_B = 2
    initial_C = 2

    # Define equations based on the methods
    # Method 1: A -> 2X + Y
    eq1 = Eq(X, 2 * A)
    eq2 = Eq(Y, A)

    # Method 2: 2B -> X
    eq3 = Eq(X, X + B // 2)

    # Method 3: 2C -> Y
    eq4 = Eq(Y, Y + C // 2)

    # Method 4: B + C -> A
    eq5 = Eq(A, A + Min(B, C))

    # Method 5: X + Y -> Z
    eq6 = Eq(Z, Min(X, Y))

    # Solve the system of equations
    solution = solve((eq1, eq2, eq3, eq4, eq5, eq6), (A, B, C, X, Y, Z), dict=True)

    # Extract the solution
    if solution:
        final_counts = solution[0]
        result = [str(final_counts[A]), str(final_counts[B]), str(final_counts[C]),
                  str(final_counts[X]), str(final_counts[Y]), str(final_counts[Z])]
    else:
        result = ["0", "0", "0", "0", "0", "0"]

    print(result)

symbolic_dismantle_and_combine()