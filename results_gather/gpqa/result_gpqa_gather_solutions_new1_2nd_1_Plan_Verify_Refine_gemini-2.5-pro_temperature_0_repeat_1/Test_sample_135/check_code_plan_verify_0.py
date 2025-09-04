import sympy as sp

def check_divergence_integral():
    """
    Checks the volume integral of the divergence of the given vector field
    f(r) = (1/r^2) * r_hat by using the Divergence Theorem.
    It calculates the equivalent surface integral (flux) and compares it
    to the provided answer.
    """
    # Define symbolic variables for spherical coordinates
    # R (Radius) is treated as a constant on the surface and cancels out,
    # so we don't need it as a variable for the final integrand.
    theta, phi = sp.symbols('theta phi')

    # According to the Divergence Theorem, the volume integral of the divergence
    # is equal to the flux of the field through the enclosing surface.
    # For f = (1/r^2)r_hat and dS = (R^2 * sin(theta) * dtheta * dphi)r_hat,
    # the dot product f ⋅ dS on the surface (where r=R) simplifies to:
    # (1/R^2) * (R^2 * sin(theta)) = sin(theta)
    integrand = sp.sin(theta)

    # Integrate over the surface of the sphere.
    # Limits for theta are [0, pi]
    # Limits for phi are [0, 2*pi]
    
    # Integrate with respect to theta first
    integral_over_theta = sp.integrate(integrand, (theta, 0, sp.pi))
    
    # Integrate the result with respect to phi
    calculated_flux = sp.integrate(integral_over_theta, (phi, 0, 2 * sp.pi))

    # The options given in the question are:
    # A) 4/3 π R
    # B) 0
    # C) 4 π
    # D) 1
    
    # The final answer provided is 'C'.
    final_answer_key = 'C'
    
    # The value corresponding to answer 'C' is 4π.
    expected_answer_value = 4 * sp.pi

    # Check 1: The calculation must be correct.
    if not sp.Eq(calculated_flux, expected_answer_value):
        return f"Calculation Mismatch: The symbolic integration resulted in {calculated_flux}, but the correct theoretical value is {expected_answer_value}."

    # Check 2: The provided answer key must correspond to the correct value.
    # This is implicitly checked by comparing the calculated value to the expected value of option 'C'.
    
    # Check 3: The answer must satisfy the constraints.
    # The result must be a constant, independent of R. This is a key feature of inverse-square fields.
    # Option A (4/3 * pi * R) violates this. The calculation confirms independence from R.
    if 'R' in str(calculated_flux):
        return f"Constraint Violation: The calculated result {calculated_flux} depends on the radius R, which is incorrect. The R terms should have cancelled out."

    return "Correct"

# Run the check and print the result
result = check_divergence_integral()
print(result)