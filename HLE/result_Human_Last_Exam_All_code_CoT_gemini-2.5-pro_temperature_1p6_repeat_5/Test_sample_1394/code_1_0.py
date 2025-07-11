def solution_printer():
    """
    Prints the general solution to the differential equation and its components.
    The solution is y^2 = C*x^2 + C^2 + C + 9, where C is an arbitrary constant.
    """
    # The final equation.
    equation = "y**2 = C*x**2 + C**2 + C + 9"

    # The numbers appearing in the equation are:
    # 2: power of y, x, and C
    # 1: implicit coefficient of C
    # 9: the constant term
    power_y = 2
    power_x = 2
    power_C = 2
    coeff_C = 1
    constant = 9

    print("The general solution is:")
    print(equation)
    print("\nEach number in the final equation is:")
    print(f"Power of y: {power_y}")
    print(f"Power of x: {power_x}")
    print(f"Highest power of C: {power_C}")
    print(f"Coefficient of C (with power 1): {coeff_C}")
    print(f"Constant term: {constant}")

solution_printer()