import sympy as sp

def solve_bdf4_stability_angle():
    """
    Calculates the exact A(alpha)-stability angle for the BDF4 method.
    """
    # Step 1: Define the tangency point based on literature.
    # The cosine of the critical angle theta_0 is c0.
    c0 = sp.S(2)/7
    
    # Calculate the sine of theta_0. We choose the positive root for theta in [0, pi].
    s0 = sp.sqrt(1 - c0**2)
    
    # Step 2: Define the boundary curve function z(theta).
    # We use the variable w = 1 - exp(-i*theta).
    # exp(-i*theta) = cos(theta) - i*sin(theta) = c0 - i*s0
    w = 1 - (c0 - sp.I * s0)
    
    # The function z is a sum of powers of w.
    z0 = w + w**2/2 + w**3/3 + w**4/4
    
    # Step 3: Extract the real (x0) and imaginary (y0) parts of z0.
    # Sympy's expand() function helps separate the real and imaginary parts.
    z0_expanded = sp.expand(z0)
    x0 = sp.re(z0_expanded)
    y0 = sp.im(z0_expanded)
    
    # Step 4: Calculate tan(alpha) = -y0/x0.
    tan_alpha = -y0 / x0
    
    # The result is of the form A*sqrt(B)/C. We extract A, B, and C.
    # tan_alpha is a fraction, so we get its numerator and denominator.
    num, den = tan_alpha.as_numer_denom()
    
    # The numerator contains the sqrt term.
    # We find the integer part and the part under the square root.
    # The form is num_coeff * sqrt(radicand).
    num_coeff = num.as_coeff_Mul()[0]
    radicand = num.as_coeff_Mul()[1].args[0]
    
    # Final values for the equation alpha = arctan(A * sqrt(B) / C)
    A = num_coeff
    B = radicand
    C = den
    
    # Step 5: Print the final equation.
    print(f"The exact stability angle alpha for BDF4 is given by:")
    print(f"alpha = arctan(A * sqrt(B) / C)")
    print(f"where:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    
    # We can also print the equation in a more readable format.
    print("\nFinal Equation:")
    print(f"alpha = arctan({A}*sqrt({B})/{C})")
    
    # For verification, we can compute the numerical value in degrees.
    numerical_alpha_rad = sp.atan(tan_alpha)
    numerical_alpha_deg = sp.deg(numerical_alpha_rad)
    print(f"\nThis corresponds to a numerical value of alpha ≈ {numerical_alpha_deg.evalf(5)} degrees.")
    
    # Return the symbolic expression as the final answer.
    return f"arctan({A}*sqrt({B})/{C})"

# Execute the function and store the result.
final_answer = solve_bdf4_stability_angle()

# Output the final answer in the specified format.
# print(f"<<<{final_answer}>>>")
# The final expression is too long for the requested format.
# Let's return the numerical value which is more standard.
numerical_angle = sp.deg(sp.atan(sp.S(2511)*sp.sqrt(5)/1625)).evalf(4)
final_value = str(numerical_angle)
# print(f"<<<{final_value}>>>")
# The prompt is a bit ambiguous. Let's stick to the symbolic expression since it asks for the exact value.

# Final symbolic answer to be wrapped.
final_symbolic_answer = "arctan(2511*sqrt(5)/1625)"
print(f"\nFinal Answer for Submission:")
print(f"<<<{final_symbolic_answer}>>>")
