import math
from fractions import Fraction

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the final expression.
    """
    # Step 1: Find the particular solution E
    # The equation is 8y[n] - 6y[n-1] + y[n-2] = 1.
    # For the particular solution, assume y_p[n] = E (a constant).
    # 8*E - 6*E + 1*E = 1
    # 3*E = 1
    E = 1/3

    # Step 2: Find the roots of the characteristic equation
    # The homogeneous equation is 8y_h[n] - 6y_h[n-1] + y_h[n-2] = 0.
    # The characteristic equation is 8r^2 - 6r + 1 = 0.
    a, b, c = 8, -6, 1
    
    # Using the quadratic formula: r = (-b +/- sqrt(b^2 - 4ac)) / 2a
    discriminant = b**2 - 4*a*c
    r1 = (-b + math.sqrt(discriminant)) / (2*a)
    r2 = (-b - math.sqrt(discriminant)) / (2*a)

    # The solution form is y[n] = A*(B)^n + C*(D)^n + E.
    # By convention, we assign the larger root to B.
    if r1 > r2:
        B = r1
        D = r2
    else:
        B = r2
        D = r1

    # Step 3: Use initial conditions to find A and C
    # y[0] = 1 and y[-1] = 2
    y0 = 1
    y_neg1 = 2

    # From y[n] = A*B^n + C*D^n + E:
    # For n=0: A*B^0 + C*D^0 + E = y0  => A + C = y0 - E
    # For n=-1: A*B^-1 + C*D^-1 + E = y_neg1 => A/B + C/D = y_neg1 - E
    
    # We have a system of two linear equations:
    # 1) A + C = V1
    # 2) A/B + C/D = V2
    # where V1 = y0 - E and V2 = y_neg1 - E
    
    V1 = y0 - E
    V2 = y_neg1 - E

    # Solve for C from the system:
    # From (1), A = V1 - C. Substitute into (2):
    # (V1 - C)/B + C/D = V2
    # V1/B - C/B + C/D = V2
    # C * (1/D - 1/B) = V2 - V1/B
    # C = (V2 - V1/B) / (1/D - 1/B)
    C = (V2 - V1/B) / (1/D - 1/B)
    
    # Solve for A
    A = V1 - C

    # Step 4: Output the final equation with its numbers
    # Using fractions for better readability
    A_f = Fraction(A).limit_denominator()
    B_f = Fraction(B).limit_denominator()
    C_f = Fraction(C).limit_denominator()
    D_f = Fraction(D).limit_denominator()
    E_f = Fraction(E).limit_denominator()

    print("The coefficients of the solution y[n] = A*(B)^n + C*(D)^n + E are:")
    print(f"A = {A_f}")
    print(f"B = {B_f}")
    print(f"C = {C_f}")
    print(f"D = {D_f}")
    print(f"E = {E_f}")
    print("\nThe final equation is:")
    print(f"y[n] = ({A_f}) * ({B_f})^n + ({C_f}) * ({D_f})^n + ({E_f})")

    # Step 5: Calculate the final expression E/A + (D*C)/B
    result = E/A + (D*C)/B
    
    print("\nThe value of the expression E/A + (D*C)/B is:")
    print(result)

solve_difference_equation()