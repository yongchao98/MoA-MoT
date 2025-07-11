from fractions import Fraction

def solve_difference_equation():
    """
    This function solves the given difference equation and calculates the required expression.
    """
    # Step 1: Find the homogeneous solution
    # The characteristic equation is 8r^2 - 6r + 1 = 0.
    # Factoring it gives (4r - 1)(2r - 1) = 0.
    # The roots are r1 = 1/2 and r2 = 1/4.
    
    # By convention, we assign the larger root to B and the smaller to D.
    B = Fraction(1, 2)
    D = Fraction(1, 4)
    
    # Step 2: Find the particular solution
    # Assume a constant particular solution y_p[n] = E.
    # Substituting into the equation: 8E - 6E + E = 1 => 3E = 1.
    E = Fraction(1, 3)
    
    # Step 3 & 4: Use initial conditions to find A and C
    # The general solution is y[n] = A*B^n + C*D^n + E
    # Initial condition y[0] = 1:
    # A*B^0 + C*D^0 + E = 1  => A + C = 1 - E
    # Initial condition y[-1] = 2:
    # A*B^-1 + C*D^-1 + E = 2 => A/B + C/D = 2 - E

    # We have a system of two linear equations:
    # 1) A + C = 1 - 1/3 = 2/3
    # 2) A/(1/2) + C/(1/4) = 2 - 1/3 => 2A + 4C = 5/3
    
    # Step 5: Solve the system
    # From (1), A = 2/3 - C. Substitute into (2):
    # 2*(2/3 - C) + 4C = 5/3
    # 4/3 - 2C + 4C = 5/3
    # 2C = 5/3 - 4/3
    # 2C = 1/3 => C = 1/6
    C = Fraction(1, 6)
    
    # Substitute C back into A = 2/3 - C
    # A = 2/3 - 1/6 = 4/6 - 1/6 = 3/6 = 1/2
    A = Fraction(1, 2)
    
    # The solution is confirmed. Let's print the values.
    print(f"The constants for the solution y[n] = A*(B)^n + C*(D)^n + E are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print("\nThe final equation is:")
    print(f"y[n] = ({A}) * ({B})^n + ({C}) * ({D})^n + ({E})")

    # Step 6: Calculate the final expression E/A + (D*C)/B
    result = E/A + (D*C)/B
    
    print("\nCalculating the expression E/A + (D*C)/B:")
    print(f"E/A = ({E}) / ({A}) = {E/A}")
    print(f"(D*C)/B = (({D}) * ({C})) / ({B}) = {(D*C)/B}")
    print(f"The final result is {E/A} + {(D*C)/B} = {result}")

if __name__ == '__main__':
    solve_difference_equation()