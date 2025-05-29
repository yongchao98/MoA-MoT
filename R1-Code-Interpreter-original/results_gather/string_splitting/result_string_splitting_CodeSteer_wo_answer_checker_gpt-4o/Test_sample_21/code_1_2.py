from sympy import symbols, Eq, solve, floor

def symbolic_dismantle_and_combine():
    # Define symbolic variables
    A, B, C, X, Y, Z = symbols('A B C X Y Z', integer=True)

    # Initial counts
    initial_A = 3
    initial_B = 2
    initial_C = 2

    # Convert B and C to A as much as possible
    convert_to_A = min(initial_B, initial_C)
    A_final = initial_A + convert_to_A
    B_final = initial_B - convert_to_A
    C_final = initial_C - convert_to_A

    # Calculate parts X and Y from A
    X_from_A = 2 * A_final
    Y_from_A = A_final

    # Calculate additional parts X and Y from remaining B and C
    X_from_B = floor(B_final / 2)
    Y_from_C = floor(C_final / 2)

    # Total parts X and Y
    total_X = X_from_A + X_from_B
    total_Y = Y_from_A + Y_from_C

    # Calculate parts Z
    total_Z = min(total_X, total_Y)

    # Remaining parts X and Y after converting to Z
    remaining_X = total_X - total_Z
    remaining_Y = total_Y - total_Z

    # Output the final counts
    result = [str(0), str(B_final % 2), str(C_final % 2), str(remaining_X), str(remaining_Y), str(total_Z)]
    print(result)

symbolic_dismantle_and_combine()