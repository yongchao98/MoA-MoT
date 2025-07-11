from fractions import Fraction

def solve_difference_equation():
    """
    Solves the difference equation 8y[n] - 6y[n-1] + y[n-2] = 1
    with initial conditions y[0] = 1 and y[-1] = 2.
    The solution form is y[n] = A*(B**n) + C*(D**n) + E.
    It then calculates the value of E/A + (D*C)/B.
    """
    print("### Solving the difference equation: 8y[n] - 6y[n-1] + y[n-2] = 1 ###")
    print("Initial conditions: y[0] = 1, y[-1] = 2")
    print("Solution form: y[n] = A*(B**n) + C*(D**n) + E\n")

    # Step 1: Find Homogeneous Solution Roots (B and D)
    # The characteristic equation is 8r^2 - 6r + 1 = 0.
    # We solve for r using the quadratic formula: r = (-b +/- sqrt(b^2 - 4ac)) / 2a
    a, b, c = 8, -6, 1
    delta = b**2 - 4*a*c
    sqrt_delta = delta**0.5

    r1 = (-b + sqrt_delta) / (2*a)
    r2 = (-b - sqrt_delta) / (2*a)

    # For consistency, assign the larger root to B and the smaller root to D.
    B = Fraction(max(r1, r2)).limit_denominator()
    D = Fraction(min(r1, r2)).limit_denominator()

    print("--- Step 1: Homogeneous Solution ---")
    print(f"The characteristic equation is 8r^2 - 6r + 1 = 0.")
    print(f"The roots are r1 = {r1} and r2 = {r2}.")
    print(f"Assigning the larger root to B and the smaller to D gives:")
    print(f"B = {B}")
    print(f"D = {D}\n")

    # Step 2: Find the Particular Solution (E)
    # For a constant input f[n]=1, we assume a constant particular solution y_p[n] = E.
    # Substituting into the original equation: 8*E - 6*E + 1*E = 1, which simplifies to 3*E = 1.
    E = Fraction(1, 3)

    print("--- Step 2: Particular Solution ---")
    print("Assuming a particular solution of the form y_p[n] = E, we get:")
    print("8*E - 6*E + E = 1  =>  3*E = 1")
    print(f"E = {E}\n")

    # Step 3: Find Coefficients A and C using initial conditions
    # The general solution is y[n] = A*(B**n) + C*(D**n) + E.
    # For n=0: y[0] = A*(B**0) + C*(D**0) + E = 1  => A + C = 1 - E
    # For n=-1: y[-1] = A*(B**-1) + C*(D**-1) + E = 2 => A/B + C/D = 2 - E
    
    # We have a system of two linear equations:
    # 1) A + C = 1 - E
    # 2) (1/B)*A + (1/D)*C = 2 - E
    val1 = 1 - E
    val2 = 2 - E
    
    # Solving the system:
    # From (1), A = val1 - C. Substitute into (2):
    # (1/B)*(val1 - C) + (1/D)*C = val2
    # C * (1/D - 1/B) = val2 - val1/B
    # C * (B - D)/(B*D) = (val2*B - val1)/B
    # C = (val2*B - val1) * D / (B - D)
    C = (val2 * B - val1) * D / (B - D)
    A = val1 - C

    print("--- Step 3: Coefficients A and C ---")
    print("Using initial conditions y[0]=1 and y[-1]=2, we solve the system:")
    print(f"1) A + C = 1 - E  => A + C = {val1}")
    print(f"2) A/B + C/D = 2 - E => ({1/B})*A + ({1/D})*C = {val2}")
    print(f"Solving the system gives:")
    print(f"A = {A}")
    print(f"C = {C}\n")

    # Step 4: State the final equation
    print("--- Step 4: The Final Closed-Form Solution ---")
    print("The full solution is y[n] = A*(B**n) + C*(D**n) + E")
    print(f"y[n] = ({A}) * ({B})^n + ({C}) * ({D})^n + ({E})\n")

    # Step 5: Calculate the required expression
    final_expression_value = E/A + (D*C)/B

    print("--- Step 5: Calculation of the Final Expression ---")
    print(f"The value of E/A is ({E}) / ({A}) = {E/A}")
    print(f"The value of (D*C)/B is (({D})*({C})) / ({B}) = {(D*C)/B}")
    print(f"The final result is E/A + (D*C)/B = {E/A} + {(D*C)/B} = {final_expression_value}")
    
    # Print the final answer as a decimal
    print(f"\nFinal Answer (as a decimal): {float(final_expression_value)}")

if __name__ == "__main__":
    solve_difference_equation()