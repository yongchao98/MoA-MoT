import scipy.integrate as integrate
import numpy as np

def calculate_fz(z_val):
    """
    Calculates the exact value of f_Z(z) at z=z_val using numerical integration.

    Args:
        z_val (float): The point at which to evaluate the PDF f_Z(z).

    Returns:
        float: The value of f_Z(z_val).
    """

    # Step 9: Define the CDF and PDF of the distance D = |X - x| for X ~ U(0,1).
    def F_x(y, x):
        """CDF F_x(y) = P(|X-x| <= y)"""
        if y < 0:
            return 0
        return max(0, min(x + y, 1) - max(x - y, 0))

    def f_x(y, x):
        """PDF f_x(y) for the distance |X-x|"""
        if y < 0:
            return 0
        if y <= min(x, 1 - x):
            return 2.0
        elif y <= max(x, 1 - x):
            return 1.0
        else:
            return 0.0

    # This is the integrand for the inner integral, which defines f_{Z|X_1}(z|x).
    # It calculates (1/y) * f_{D_(2)|X_1=x}(y).
    def inner_integrand(y, x, z):
        """The function to integrate over y."""
        if y == 0:
            return 0
        
        Fx_y = F_x(y, x)
        fx_y = f_x(y, x)
        
        # PDF of the second order statistic of 3 distances
        f_D2_y_given_x = 6.0 * Fx_y * (1.0 - Fx_y) * fx_y
        
        return (1.0 / y) * f_D2_y_given_x

    # Step 8: Define the conditional PDF f_{Z|X_1}(z|x) as an integral.
    def conditional_pdf_f_z_given_x(x, z):
        """Computes f_{Z|X_1=x}(z) via numerical integration."""
        lower_bound_y = abs(z - x)
        upper_bound_y = max(x, 1 - x)
        
        if lower_bound_y >= upper_bound_y:
            return 0.0

        integral_val, _ = integrate.quad(
            inner_integrand, lower_bound_y, upper_bound_y, args=(x, z)
        )
        return integral_val

    # Step 5: Calculate f_Z(z) by integrating f_{Z|X_1}(z|x) over x from 0 to 1.
    final_value, _ = integrate.quad(
        conditional_pdf_f_z_given_x, 0, 1, args=(z_val,)
    )
    
    return final_value

# Set the value of z for which we want to calculate the PDF.
z = 0.2
result = calculate_fz(z)

# The result is the exact value requested.
print(f"To calculate f_Z({z}), we solve the double integral:")
print(f"f_Z({z}) = integral from x=0 to 1 of [ integral from y=|{z}-x| to max(x,1-x) of g(y,x) dy ] dx")
print("where g(y,x) is derived from the PDF of the second closest distance.")
print(f"The exact value is: {result}")
print("\nThe final equation with the calculated value is:")
print(f"f({z}) = {result}")
