from sympy import symbols, Eq, solve, floor

def symbolic_dismantle_and_combine():
    # Define symbolic variables
    A, B, C, X, Y, Z = symbols('A B C X Y Z', integer=True)

    # Initial counts
    initial_A = 3
    initial_B = 2
    initial_C = 2

    # Initialize the counts
    machines_A = initial_A
    machines_B = initial_B
    machines_C = initial_C
    parts_X = 0
    parts_Y = 0
    parts_Z = 0

    # Flag to check if any method was executed in the last cycle
    any_method_executed = True

    while any_method_executed:
        any_method_executed = False

        # Method 1: A -> 2X + Y
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
            any_method_executed = True

        # Method 2: 2B -> X
        if machines_B >= 2:
            machines_B -= 2
            parts_X += 1
            any_method_executed = True

        # Method 3: 2C -> Y
        if machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
            any_method_executed = True

        # Method 4: B + C -> A
        if machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
            any_method_executed = True

        # Method 5: X + Y -> Z
        if parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
            any_method_executed = True

    # Output the final counts
    result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
    print(result)

symbolic_dismantle_and_combine()