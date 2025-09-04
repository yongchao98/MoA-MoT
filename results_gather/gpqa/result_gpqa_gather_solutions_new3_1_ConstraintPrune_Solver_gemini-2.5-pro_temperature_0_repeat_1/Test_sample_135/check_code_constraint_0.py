import sympy

def check_integral_correctness():
    """
    Verifies the solution to the volume integral of the divergence problem
    by symbolically calculating the equivalent surface integral via the Divergence Theorem.

    The problem is to evaluate ∫_V (∇ ⋅ f) dV for f(r) = (1/r²) * r̂.
    The Divergence Theorem transforms this into ∮_S (f ⋅ dS).

    On the surface of a sphere of radius R:
    - f = (1/R²) * r̂
    - dS = (R² * sin(θ) * dθ * dφ) * r̂
    - The dot product f ⋅ dS simplifies to sin(θ) dθ dφ.

    The code calculates ∫[φ=0 to 2π] ∫[θ=0 to π] sin(θ) dθ dφ.
    """
    
    # 1. Define the symbolic variables for the angles in spherical coordinates.
    theta, phi = sympy.symbols('theta phi')

    # 2. Define the integrand, which is sin(θ) after simplification.
    integrand = sympy.sin(theta)

    # 3. Perform the definite double integral over the surface of a sphere.
    # The limits for the polar angle θ are from 0 to π.
    # The limits for the azimuthal angle φ are from 0 to 2π.
    try:
        # Integrate with respect to theta first.
        integral_over_theta = sympy.integrate(integrand, (theta, 0, sympy.pi))
        
        # Then, integrate the result with respect to phi.
        calculated_value = sympy.integrate(integral_over_theta, (phi, 0, 2 * sympy.pi))
    except Exception as e:
        return f"An error occurred during the symbolic integration: {e}"

    # 4. The final answer provided in the prompt is 'A', which corresponds to 4π.
    # Let's define the expected value based on the provided answer.
    expected_value = 4 * sympy.pi

    # 5. Check if the calculated value matches the expected value.
    # We use sympy.simplify to ensure a robust comparison.
    if sympy.simplify(calculated_value - expected_value) == 0:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer 'A' corresponds to the value {expected_value}. "
                f"However, the symbolic calculation of the integral ∮sin(θ)dθdφ yields {calculated_value}. "
                f"The provided answer does not match the calculated result.")

# Execute the check
result = check_integral_correctness()
print(result)