import sympy
from sympy import sin, pi, integrate, symbols

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer by re-evaluating the integral using symbolic mathematics.

    The problem is to evaluate the volume integral of the divergence of a vector field:
    ∫_V (∇ ⋅ f) dV
    where V is a sphere of radius R centered at the origin, and the field is f(r) = (1/r^2) * r_hat.

    The provided solution correctly applies the Divergence Theorem to convert the volume integral
    into a surface integral over the boundary of the sphere (S):
    ∫_V (∇ ⋅ f) dV = ∮_S (f ⋅ dS)

    This function will calculate the surface integral ∮_S (f ⋅ dS) symbolically to verify the result.
    """
    try:
        # 1. Define symbolic variables for the integration.
        # We use the angular variables theta (polar angle) and phi (azimuthal angle).
        theta, phi = symbols('theta phi', real=True)

        # 2. Set up the integrand for the surface integral.
        # On the surface of the sphere, the radius r is a constant R.
        # The vector field f is f(R) = (1/R^2) * r_hat.
        # The differential surface area element in spherical coordinates is dS = R^2 * sin(theta) * dtheta * dphi * r_hat.
        # The direction of dS is radially outward, parallel to r_hat.

        # The dot product f ⋅ dS is:
        # f ⋅ dS = ( (1/R^2) * r_hat ) ⋅ ( R^2 * sin(theta) * dtheta * dphi * r_hat )
        # Since r_hat ⋅ r_hat = 1, and the R^2 terms cancel out, the expression to integrate is:
        # sin(theta) * dtheta * dphi
        # So, the integrand for our double integral is sin(theta).
        integrand = sin(theta)

        # 3. Perform the double integration over the surface of the sphere.
        # The polar angle theta ranges from 0 to pi.
        # The azimuthal angle phi ranges from 0 to 2*pi.

        # First, integrate with respect to theta from 0 to pi.
        integral_over_theta = integrate(integrand, (theta, 0, pi))

        # Then, integrate the result with respect to phi from 0 to 2*pi.
        surface_integral_result = integrate(integral_over_theta, (phi, 0, 2*pi))

        # 4. Compare the calculated result with the given answer.
        # The answer from the LLM is C, which corresponds to the value 4π.
        expected_answer_value = 4 * pi

        # Use sympy.equals() for a robust symbolic comparison.
        if sympy.equals(surface_integral_result, expected_answer_value):
            return "Correct"
        else:
            return (f"The answer is incorrect. "
                    f"The calculation of the surface integral using the Divergence Theorem yields {surface_integral_result}, "
                    f"but the provided answer corresponds to {expected_answer_value}.")

    except Exception as e:
        return f"An error occurred during the verification process: {e}"

# Run the check
result = check_answer_correctness()
print(result)