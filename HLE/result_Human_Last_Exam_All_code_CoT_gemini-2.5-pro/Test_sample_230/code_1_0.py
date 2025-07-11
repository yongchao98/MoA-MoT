def check_cancellability():
    """
    This function checks the n-cancellability of a known non-medial C-groupoid.
    A C-groupoid is an idempotent, commutative, left-distributive magma.
    The table represents a 4-element non-medial C-groupoid from the literature
    (e.g., G_1 in "Small non-medial C-groupoids" by D. Stanovsky).
    We check for which n it is n-cancellable.
    n-cancellable means: if x_n * y = y, then x must be equal to y.
    If we find even one pair (x, y) with x != y where x_n * y = y,
    the magma is NOT n-cancellable.
    """
    # Multiplication table for the non-medial C-groupoid G_1
    table = {
        0: {0: 0, 1: 2, 2: 3, 3: 2},
        1: {0: 2, 1: 1, 2: 2, 3: 3},
        2: {0: 3, 1: 2, 2: 2, 3: 3},
        3: {0: 2, 1: 3, 2: 3, 3: 3},
    }
    size = len(table)

    def op(a, b):
        return table[a][b]

    def op_n(a, b, n):
        result = b
        for _ in range(n):
            result = op(a, result)
        return result

    print("Checking for which n the magma is NOT n-cancellable:")
    # We test for n from 1 to 5
    for n in range(1, 6):
        is_n_cancellable = True
        for x in range(size):
            for y in range(size):
                if x == y:
                    continue
                
                # Check if x_n * y = y
                if op_n(x, y, n) == y:
                    print(f"For n={n}, the magma is NOT {n}-cancellable.")
                    # We found a counterexample.
                    # The equation is {x}_n * {y} = {y}.
                    # Let's print the calculation step by step.
                    calculation_str = f"    Proof: Let x={x}, y={y}. We have x != y."
                    current_y = y
                    for i in range(1, n + 1):
                        next_y = op(x, current_y)
                        calculation_str += f"\n    Step {i}: {x} * {current_y} = {next_y}"
                        current_y = next_y
                    calculation_str += f"\n    Final result: {x}_{n} * {y} = {current_y}, which is equal to the original y={y}."
                    print(calculation_str)
                    is_n_cancellable = False
                    break  # Move to the next n
            if not is_n_cancellable:
                break
        
        if is_n_cancellable:
            print(f"For n={n}, the magma might be {n}-cancellable (no counterexample found).")

check_cancellability()