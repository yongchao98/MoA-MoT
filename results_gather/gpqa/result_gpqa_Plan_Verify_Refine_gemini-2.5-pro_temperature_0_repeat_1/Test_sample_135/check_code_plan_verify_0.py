import sympy
from sympy import sin, pi, integrate

def check_correctness():
    """
    This function checks the correctness of the provided answer by re-calculating the integral
    using the Divergence Theorem, as described in the solution.

    The problem is to evaluate ∫_V (∇ ⋅ f) dV for f(r) = (1/r^2) * r_hat.
    The Divergence Theorem states: ∫_V (∇ ⋅ f) dV = ∮_S (f ⋅ dS).

    The check will verify the calculation of the surface integral ∮_S (f ⋅ dS) over a sphere
    of radius R.
    """
    # Define symbolic variables for the spherical coordinate angles
    theta, phi = sympy.symbols('theta phi')
    
    # According to the provided solution, the surface integral is set up as follows:
    # The vector field on the surface of the sphere (where r=R) is f = (1/R^2) * r_hat.
    # The differential surface area element is dS = R^2 * sin(θ) * dθ * dφ * r_hat.
    # The dot product f ⋅ dS simplifies to:
    # (1/R^2) * r_hat ⋅ (R^2 * sin(θ) * dθ * dφ * r_hat) = sin(θ) * dθ * dφ.
    # The integrand for the surface integral is therefore sin(θ).
    
    integrand = sin(theta)
    
    # The limits of integration for a full sphere are:
    # θ (polar angle) from 0 to π
    # φ (azimuthal angle) from 0 to 2π
    
    # The expected answer from the provided solution is 4π.
    expected_answer = 4 * pi
    
    try:
        # We perform the double integral: ∫_0^(2π) dφ ∫_0^π sin(θ) dθ
        
        # First, calculate the inner integral with respect to θ
        inner_integral = integrate(integrand, (theta, 0, pi))
        
        # The result of the inner integral should be 2. Let's verify this step.
        # ∫_0^π sin(θ) dθ = [-cos(θ)]_0^π = -cos(π) - (-cos(0)) = -(-1) - (-1) = 2
        if inner_integral != 2:
            return f"Incorrect intermediate step: The integral of sin(theta) from 0 to pi should be 2, but the symbolic calculation yielded {inner_integral}."
            
        # Now, calculate the outer integral with respect to φ
        final_result = integrate(inner_integral, (phi, 0, 2 * pi))
        
        # Compare the symbolic result with the expected answer
        if final_result == expected_answer:
            # The calculation is correct. The answer D) 4π is supported by this calculation.
            # The solution correctly applies the Divergence Theorem (despite the singularity at r=0,
            # which is a known special case where this method still yields the correct physical answer,
            # equivalent to using a Dirac delta function for the divergence).
            return "Correct"
        else:
            return f"The final calculated value of the integral is {final_result}, which does not match the provided answer's result of {expected_answer}. The provided answer's calculation is correct, but this check failed."

    except Exception as e:
        return f"An error occurred during the symbolic calculation: {e}"

# Run the check and print the result
result = check_correctness()
print(result)