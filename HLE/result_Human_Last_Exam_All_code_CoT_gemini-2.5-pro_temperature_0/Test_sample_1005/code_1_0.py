import sympy

def solve_integral():
    """
    This function demonstrates the analytical solution of the definite integral
    by using symbolic mathematics to show the cancellation of complex terms.
    """
    # Define the symbolic variable
    theta = sympy.Symbol('theta')
    
    print("Step 1: The integral is split into I1 and I2 and transformed using trigonometric substitutions.")
    print("The new integration variable is theta, with limits from 0 to pi/4 for both parts.\n")

    # After substitution, I1 becomes C1 * Integral(sin(theta)**(1/4), dtheta)
    C1 = 2**sympy.Rational(-15, 16)
    
    # After substitution, I2 becomes C2 * Integral(sin(theta)**(1/4) * sec(theta)**2, dtheta)
    C2 = 2**sympy.Rational(17, 16)

    # Define the common unevaluated integral part
    common_integral = sympy.Integral(sympy.sin(theta)**sympy.Rational(1, 4), (theta, 0, sympy.pi/4))

    print("The first part of the integral, I1, transforms to:")
    I1_expr = C1 * common_integral
    sympy.pprint(I1_expr)
    print("-" * 30)

    print("Step 2: We use integration by parts for the second integral, I2.")
    print("Let u = sin(theta)^(1/4) and dv = sec(theta)^2 d(theta).")
    print("The formula is Integral(u*dv) = [u*v] - Integral(v*du).")
    print("This shows that Integral(sin(theta)^(1/4)*sec(theta)^2) = [sin(theta)^(1/4)*tan(theta)] - 1/4 * Integral(sin(theta)^(1/4)).\n")

    # Evaluate the [u*v] part from 0 to pi/4
    # u*v at pi/4 is (sin(pi/4))^(1/4) * tan(pi/4) = (1/sqrt(2))^(1/4) * 1 = 2**(-1/8)
    # u*v at 0 is 0.
    evaluated_uv_simplified = 2**sympy.Rational(-1, 8)

    # Express I2 using the result of integration by parts
    I2_expr = C2 * (evaluated_uv_simplified - sympy.Rational(1, 4) * common_integral)
    
    print("The second part of the integral, I2, becomes:")
    sympy.pprint(I2_expr)
    print("-" * 30)

    print("Step 3: The total integral I = I1 + I2 is the sum of the two expressions.")
    I_total_expr = I1_expr + I2_expr
    
    # Expand the expression to show the cancellation
    I_total_expanded = sympy.expand(I_total_expr)
    
    print("The full expression for I is:")
    sympy.pprint(I_total_expanded)
    print("\nAs we can see, the integral terms have coefficients that are equal and opposite, so they cancel out.")
    print(f"Coefficient of integral in I1: {C1}")
    print(f"Coefficient of integral in I2: {-(C2/4).simplify()}")
    print("-" * 30)

    print("Step 4: The final result is the remaining term.")
    final_result = C2 * evaluated_uv_simplified
    
    # The final equation is I = 2^(17/16) * 2^(-1/8)
    print("The final equation is:")
    print(f"I = ({C2}) * ({evaluated_uv_simplified})")
    
    # Simplify the expression
    final_symbolic_answer = final_result.simplify()
    
    print("\nSimplifying the expression gives:")
    final_equation_str = f"I = {final_symbolic_answer}"
    print(final_equation_str)
    
    # Output the numbers in the final equation
    base, exp = final_symbolic_answer.as_base_exp()
    num, den = exp.as_numer_denom()
    print("\nNumbers in the final equation:")
    print(f"Base: {base}")
    print(f"Exponent Numerator: {num}")
    print(f"Exponent Denominator: {den}")

    # Calculate the numerical value
    final_numerical_value = float(final_symbolic_answer)
    print(f"\nThe numerical value of the integral is: {final_numerical_value}")
    
    return final_numerical_value

if __name__ == '__main__':
    solve_integral()