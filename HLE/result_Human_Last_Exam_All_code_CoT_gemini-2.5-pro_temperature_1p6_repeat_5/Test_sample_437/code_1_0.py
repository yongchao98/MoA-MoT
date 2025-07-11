import sympy as sp

def evaluate_1s_slater_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a Slater-type orbital.
    """
    # Define the symbols used in the calculation.
    # r, theta, phi are the spherical coordinates.
    # zeta is the exponent of the Slater orbital.
    r, theta, phi = sp.symbols('r theta phi', real=True, positive=True)
    zeta = sp.symbols('zeta', real=True, positive=True, nonzero=True)
    N = sp.Symbol('N', real=True, positive=True) # Normalization constant

    print("Step 1: Define the unnormalized 1s Slater orbital.")
    # The functional form of a 1s STO is exp(-zeta*r).
    phi_unnormalized = sp.exp(-zeta * r)
    print(f"phi_unnormalized = {phi_unnormalized}\n")

    print("Step 2: Calculate the normalization constant N.")
    print("The normalization condition is integral( (N*phi)^2 * d_tau ) = 1.")
    # The integrand for normalization includes the square of the orbital and the volume element.
    # d_tau = r^2 * sin(theta) * dr * dtheta * dphi
    integrand_norm = phi_unnormalized**2 * r**2 * sp.sin(theta)
    
    # Perform the integration over all space.
    # The integral over phi (0 to 2*pi) of sin(theta) gives 2*pi.
    # The integral over theta (0 to pi) of sin(theta) gives 2.
    # Total angular integral is 4*pi.
    # The radial integral is integral(r^2 * exp(-2*zeta*r), r from 0 to oo)
    # The total integral is N^2 * (4*pi) * integral_r = 1
    
    # Using sympy to do it step-by-step:
    norm_integral_val = sp.integrate(integrand_norm, (phi, 0, 2*sp.pi), (theta, 0, sp.pi), (r, 0, sp.oo))
    
    print(f"The integral of (phi_unnormalized)^2 over all space is: {norm_integral_val}")
    
    # Solve for N from N**2 * norm_integral_val = 1
    N_squared = 1 / norm_integral_val
    N_expr = sp.sqrt(N_squared)
    print(f"The normalization constant N is sqrt(1 / ({norm_integral_val})) = {N_expr}\n")
    
    # Define the fully normalized 1s Slater orbital
    phi_1s = N_expr * phi_unnormalized

    print("Step 3: Set up the integral for the expectation value <phi_1s| 1/r |phi_1s>.")
    operator = 1 / r
    
    # The integrand for the expectation value.
    integrand_exp_val = phi_1s * operator * phi_1s * r**2 * sp.sin(theta)
    
    print(f"Normalized phi_1s = {phi_1s}")
    print(f"Operator = {operator}")
    print("Integrand = (phi_1s)^2 * (1/r) * r^2 * sin(theta)\n")

    print("Step 4: Evaluate the integral.")
    # Perform the full integration.
    final_result = sp.integrate(integrand_exp_val, (phi, 0, 2*sp.pi), (theta, 0, sp.pi), (r, 0, sp.oo))
    
    print("The integral <phi_1s| 1/r |phi_1s> evaluates to:")
    # Print the final symbolic result as an equation
    final_eq = sp.Eq(sp.Symbol('<phi_1s|1/r|phi_1s>'), final_result)
    print(final_eq)
    print("\n----------------------------------------------------")
    print("Final Equation Analysis:")
    # Decompose the result into its numerical coefficient and variable parts
    # This addresses the prompt to "output each number in the final equation"
    coeff, terms = final_result.as_coeff_Mul()
    variable_part = 1
    if terms:
        variable_part = sp.Mul(*terms)

    print("The final result can be written as: Coefficient * Variable(s)")
    print(f"The numerical coefficient in the equation is: {coeff}")
    print(f"The variable part of the equation is: {variable_part}")

evaluate_1s_slater_integral()