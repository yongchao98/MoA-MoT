def solve_hitting_time_moment(c: float):
    """
    Calculates the supremum of alpha for which the alpha-th moment of the hitting time is finite.

    Args:
        c: A positive constant from the problem description.
    """
    if c <= 0:
        print("The constant c must be positive.")
        return

    # The formula for the supremum of alpha is (1 + 4*c) / 2
    # The numbers in the equation are:
    val_1 = 1
    val_4 = 4
    val_2 = 2

    # Calculate the supremum of alpha
    sup_alpha = (val_1 + val_4 * c) / val_2

    # Print the equation with all its numbers and the final result
    print(f"The problem is to find sup{{alpha: E[tau^alpha] < infinity}}.")
    print(f"The supremum of alpha is given by the formula (1 + 4*c) / 2.")
    print(f"For the given value c = {c}, the calculation is:")
    print(f"({val_1} + {val_4} * {c}) / {val_2} = {sup_alpha}")

# --- User Input ---
# You can change the value of c here.
# For example, let's use a sample value of c = 1.25
c_value = 1.25

solve_hitting_time_moment(c_value)
