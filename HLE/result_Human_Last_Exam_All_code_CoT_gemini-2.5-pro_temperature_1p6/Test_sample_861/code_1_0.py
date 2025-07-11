import math
from fractions import Fraction

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the final expression.
    """
    # The difference equation is 8y[n] - 6y[n-1] + y[n-2] = 1
    # The general solution is y[n] = A*(B**n) + C*(D**n) + E
    # Initial conditions are y[0] = 1, y[-1] = 2

    # Step 1: Find the particular solution, E
    # For a constant input, we assume a constant particular solution y_p[n] = E.
    # Substitute into the equation: 8*E - 6*E + E = 1 => 3*E = 1
    E = Fraction(1, 3)

    # Step 2: Find the roots of the characteristic equation for the homogeneous part
    # 8*r^2 - 6*r + 1 = 0
    a, b, c = 8, -6, 1
    # Using the quadratic formula: r = (-b +/- sqrt(b^2 - 4ac)) / 2a
    discriminant = b**2 - 4 * a * c
    sqrt_discriminant = math.sqrt(discriminant)
    
    r1 = (-b + sqrt_discriminant) / (2 * a)
    r2 = (-b - sqrt_discriminant) / (2 * a)

    # Convert roots to Fractions for precision
    root1 = Fraction(r1)
    root2 = Fraction(r2)

    # By convention, assign the larger root to B and the smaller to D
    if root1 > root2:
        B = root1
        D = root2
    else:
        B = root2
        D = root1

    # Step 3: Use initial conditions to find A and C
    # The total solution is y[n] = A*(B**n) + C*(D**n) + E
    # We have two initial conditions: y[0] = 1 and y[-1] = 2
    y0 = Fraction(1)
    y_minus_1 = Fraction(2)

    # Equation for n=0: y[0] = A*B**0 + C*D**0 + E => y[0] = A + C + E
    # A + C = y[0] - E
    eq1_rhs = y0 - E

    # Equation for n=-1: y[-1] = A*B**(-1) + C*D**(-1) + E => y[-1] = A/B + C/D + E
    # A/B + C/D = y[-1] - E
    eq2_rhs = y_minus_1 - E

    # We have a system of 2 linear equations for A and C:
    # 1) A + C = eq1_rhs
    # 2) (1/B)*A + (1/D)*C = eq2_rhs
    
    # Let's solve the system. From (1), C = eq1_rhs - A.
    # Substitute into (2): (1/B)*A + (1/D)*(eq1_rhs - A) = eq2_rhs
    # A * (1/B - 1/D) = eq2_rhs - eq1_rhs / D
    # A * (D - B)/(B*D) = (eq2_rhs*D - eq1_rhs)/D
    # A = ((eq2_rhs*D - eq1_rhs)/D) * (B*D)/(D - B) = (eq2_rhs*D - eq1_rhs) * B/(D - B)
    
    # We found A=1/2 and C=1/6 from our scratchpad calculation. Let's verify.
    A_numerator = eq2_rhs * D - eq1_rhs
    A_denominator = D - B
    A = (A_numerator / A_denominator) * B

    C = eq1_rhs - A

    # Step 4: Print the final equation with all its numbers
    print("The solved constants for the equation y[n] = A*(B**n) + C*(D**n) + E are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print("\nThus, the final equation is:")
    print(f"y[n] = ({A}) * ({B})^n + ({C}) * ({D})^n + ({E})")

    # Step 5: Calculate the final required expression
    # Expression: E/A + (D*C)/B
    term1 = E / A
    term2 = (D * C) / B
    result = term1 + term2
    
    print("\nCalculating the expression E/A + (D*C)/B:")
    print(f"The value of the expression is {result}")
    # You can also display it as a float
    # print(f"As a decimal, the value is {float(result)}")

solve_difference_equation()