import sympy

def solve_spaceship_distance():
    """
    Calculates the symbolic formula for the maximum distance l_max
    of a spaceship from an asteroid using the work-energy theorem.
    """
    # Define the symbolic variables
    m, M, l_0, v_0, F, G, l_max = sympy.symbols('m M l_0 v_0 F G l_max', real=True, positive=True)

    # From the work-energy theorem, we derived the quadratic equation for l_max:
    # F*l_max**2 - (F*l_0 + G*m*M/l_0 - 1/2*M*v_0**2)*l_max + G*m*M = 0
    # This is in the form a*x**2 + b*x + c = 0

    # Define the coefficients a, b, and c
    a = F
    b = -(F * l_0 + (G * m * M) / l_0 - sympy.Rational(1, 2) * M * v_0**2)
    c = G * m * M

    # Use the quadratic formula to find the solutions for l_max.
    # The two possible solutions are (-b +/- sqrt(b**2 - 4*a*c)) / (2*a)
    # The maximum distance corresponds to the larger root (the '+' case).
    
    discriminant = b**2 - 4 * a * c
    
    solution_l_max = (-b + sympy.sqrt(discriminant)) / (2 * a)

    # To make the output more readable, let's define an intermediate term C
    # where C = -b
    C = F * l_0 + (G * m * M) / l_0 - sympy.Rational(1, 2) * M * v_0**2
    
    # Print the result in a readable format
    print("The maximum distance l_max is given by the formula:")
    print("l_max = (C + sqrt(C**2 - 4*F*G*m*M)) / (2*F)")
    print("\nwhere the terms in the equation are:")
    print(f"C = {C}")
    print(f"F = {F}")
    print(f"G = {G}")
    print(f"m = {m}")
    print(f"M = {M}")
    
    # We can also print the full expanded expression
    # print("\nFully expanded expression for l_max:")
    # sympy.pprint(solution_l_max, use_unicode=False)


if __name__ == '__main__':
    solve_spaceship_distance()
