import math
from fractions import Fraction

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the final expression.
    """
    print("Step 1 & 2: Finding roots (B, D) and the particular solution (E)")
    # Homogeneous equation: 8y[n] - 6y[n-1] + y[n-2] = 0
    # Characteristic equation: 8r^2 - 6r + 1 = 0
    a, b, c = 8, -6, 1
    
    # Solve quadratic equation for roots
    discriminant = b**2 - 4*a*c
    r1 = (-b + math.sqrt(discriminant)) / (2*a)
    r2 = (-b - math.sqrt(discriminant)) / (2*a)
    
    # Assign roots to B and D using Fraction for precision. Let B be the larger root.
    B = Fraction(r1).limit_denominator()
    D = Fraction(r2).limit_denominator()
    
    # Particular solution for y_p[n] = E
    # 8E - 6E + E = 1 => 3E = 1
    E = Fraction(1, 3)
    
    print(f"The roots of the characteristic equation are {B} and {D}.")
    print(f"Let B = {B} and D = {D}.")
    print(f"The particular solution constant E = {E}.")
    print("-" * 20)

    print("Step 3 & 4: Solving for coefficients A and C")
    # General solution: y[n] = A*(B)^n + C*(D)^n + E
    # Initial conditions: y[0] = 1, y[-1] = 2
    
    # Using y[0] = 1:
    # A*(B)^0 + C*(D)^0 + E = 1  => A + C + E = 1 => A + C = 1 - E
    # Using y[-1] = 2:
    # A*(B)^-1 + C*(D)^-1 + E = 2 => A/B + C/D = 2 - E

    y0 = Fraction(1)
    y_minus_1 = Fraction(2)
    
    # System of equations:
    # 1) A + C = y0 - E
    # 2) A/B + C/D = y_minus_1 - E
    
    val1 = y0 - E
    val2 = y_minus_1 - E
    
    # We solve the system:
    # A + C = val1
    # (1/B)A + (1/D)C = val2
    # From eq1, A = val1 - C. Substitute into eq2:
    # (1/B)(val1 - C) + (1/D)C = val2
    # val1/B - C/B + C/D = val2
    # C * (1/D - 1/B) = val2 - val1/B
    # C = (val2 - val1/B) / (1/D - 1/B)

    C = (val2 - val1 / B) / (1/D - 1/B)
    A = val1 - C
    
    print("Using initial conditions y[0]=1 and y[-1]=2, we solve for A and C.")
    print(f"The determined coefficients are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print("-" * 20)
    
    print("The final closed-form solution is:")
    print(f"y[n] = ({A}) * ({B})^n + ({C}) * ({D})^n + ({E})")
    print("-" * 20)
    
    print("Step 5: Calculating the final expression E/A + (D*C)/B")
    
    # Calculate the expression
    result = E/A + (D*C)/B
    
    print(f"E/A = ({E}) / ({A}) = {E/A}")
    print(f"(D*C)/B = (({D}) * ({C})) / ({B}) = {(D*C)/B}")
    print(f"E/A + (D*C)/B = {E/A} + {(D*C)/B} = {result}")
    
    # Final numerical value
    final_value = float(result)
    print(f"\nThe final calculated value is {result}, which is {final_value} in decimal.")
    
    return final_value

if __name__ == '__main__':
    final_answer = solve_difference_equation()
    # The final output format expected
    # print(f"<<<{final_answer}>>>")

solve_difference_equation()
<<<0.75>>>