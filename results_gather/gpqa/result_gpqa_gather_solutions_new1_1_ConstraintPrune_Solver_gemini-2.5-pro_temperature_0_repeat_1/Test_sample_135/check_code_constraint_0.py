import sympy
from sympy import sin, pi, integrate, symbols

def check_answer_correctness():
    """
    This function verifies the answer to the physics problem by applying the Divergence Theorem.

    Problem: Evaluate the volume integral of the divergence of the field f(r) = (1/r^2) * r_hat
    inside a sphere of radius R centered at the origin.

    Method: The Divergence Theorem states that the volume integral of the divergence of a vector field
    is equal to the flux of the field through the enclosing surface.
    ∫∫∫_V (∇ ⋅ f) dV = ∬_S (f ⋅ n) dS

    The code calculates the surface integral (flux) analytically.
    """

    # 1. Define symbolic variables for the integration.
    # R is the radius of the sphere, theta is the polar angle, phi is the azimuthal angle.
    R, theta, phi = symbols('R theta phi', real=True, positive=True)

    # 2. Define the integrand for the surface integral (flux).
    # On the surface of the sphere (where r=R), the field f has magnitude 1/R^2.
    # The outward normal vector 'n' is the radial unit vector.
    # The dot product (f ⋅ n) is therefore (1/R^2).
    # The scalar surface area element dS is R^2 * sin(theta) * dtheta * dphi.
    # The full term to integrate is (f ⋅ n) * dS_scalar.
    integrand = (1 / R**2) * (R**2 * sin(theta))
    
    # The R**2 terms cancel out, simplifying the integrand.
    simplified_integrand = sympy.simplify(integrand)
    
    # At this point, simplified_integrand is sin(theta).

    # 3. Perform the double integration over the surface of the sphere.
    # The polar angle theta ranges from 0 to pi.
    # The azimuthal angle phi ranges from 0 to 2*pi.
    integral_over_theta = integrate(simplified_integrand, (theta, 0, pi))
    calculated_result = integrate(integral_over_theta, (phi, 0, 2*pi))

    # 4. Define the options and the provided answer from the prompt.
    # The question's options are: A) 1, B) 4 π, C) 0, D) 4/3 π R
    # The final consolidated answer given in the prompt is <<<B>>>.
    
    # The value corresponding to option B is 4 * pi.
    expected_value_from_answer = 4 * pi

    # 5. Compare the calculated result with the value from the provided answer.
    # The .equals() method is used for robust comparison of symbolic expressions.
    if calculated_result.equals(expected_value_from_answer):
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is 'B', which corresponds to the value {expected_value_from_answer}. "
                f"However, the correct value of the integral, calculated using the Divergence Theorem, is {calculated_result}.")

# Execute the check and print the result.
result = check_answer_correctness()
print(result)