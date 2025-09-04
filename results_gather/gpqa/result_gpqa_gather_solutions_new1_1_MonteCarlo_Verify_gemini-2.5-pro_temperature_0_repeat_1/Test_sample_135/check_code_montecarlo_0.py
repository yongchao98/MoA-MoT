import numpy as np
from scipy import integrate

def check_integral_correctness():
    """
    Verifies the solution by numerically calculating the surface integral
    derived from the Divergence Theorem.
    """
    # The problem is to evaluate the volume integral of the divergence of f = (1/r^2) * r_hat.
    # By the Divergence Theorem, ∫_V (∇ ⋅ f) dV = ∮_S (f ⋅ dS).
    # This surface integral simplifies to ∫∫ sin(θ) dθ dφ over the surface of a sphere.
    # The limits are θ from 0 to π, and φ from 0 to 2π.

    # We will use scipy's double integrator (dblquad) to compute this.
    # The function to integrate is sin(θ). dblquad expects a function of two
    # variables, f(y, x), where y is the variable for the inner integral.
    # Let's set y = θ and x = φ.
    integrand = lambda theta, phi: np.sin(theta)

    # Define the integration limits for the outer integral (phi).
    phi_min = 0
    phi_max = 2 * np.pi

    # Define the integration limits for the inner integral (theta).
    # These are constant, so we can use lambda functions that return constants.
    theta_min = lambda phi: 0
    theta_max = lambda phi: np.pi

    # Perform the double integration.
    # The order of arguments is: func, x_min, x_max, y_min, y_max
    # This corresponds to ∫[x_min to x_max] ∫[y_min(x) to y_max(x)] f(y, x) dy dx
    numerical_result, error = integrate.dblquad(integrand, phi_min, phi_max, theta_min, theta_max)

    # The analytical result is known to be 4π.
    analytical_result = 4 * np.pi

    # The question provides options, and the proposed answer is 'C', which corresponds to 4π.
    # Let's check if our numerical calculation confirms this.
    if not np.isclose(numerical_result, analytical_result):
        return (f"Incorrect. The numerical calculation gives {numerical_result:.6f}, "
                f"which does not match the expected analytical result of 4π ({analytical_result:.6f}). "
                f"There might be an error in the analytical derivation or the numerical setup.")

    # The provided answer is 'C', which corresponds to the value 4π.
    # Since our calculation confirms the value is 4π, the answer 'C' is correct.
    return "Correct"

# Run the check
result_message = check_integral_correctness()
print(result_message)