import cmath
import math

def solve_bdf4_angle():
    """
    Calculates the A(alpha)-stability angle for the BDF4 method.
    """
    # The critical angle theta is pi/3, so xi = e^(i*pi/3)
    theta = math.pi / 3
    xi = cmath.exp(1j * theta)

    # BDF4 characteristic polynomial rho(xi)
    # rho(xi) = (25/12)xi^4 - 4xi^3 + 3xi^2 - (4/3)xi + 1/4
    rho_xi = ( (25/12) * xi**4 - 4 * xi**3 + 3 * xi**2 - (4/3) * xi + 1/4 )

    # The stability boundary is z(xi) = rho(xi) / sigma(xi) where sigma(xi)=xi^4
    z_t = rho_xi / (xi**4)

    # The angle alpha is |arg(-z_t)|
    # We want to find the value T such that alpha = arctan(T)
    # arg(Z) = arctan(Im(Z)/Re(Z)).
    # So alpha = |arg(-z_t)| = arctan(|-Im(z_t)/-Re(z_t)|) = arctan(|Im(z_t)/Re(z_t)|)
    re_z = z_t.real
    im_z = z_t.imag

    tangent_value = abs(im_z / re_z)
    
    # Calculate the final angle alpha in radians
    alpha_rad = math.atan(tangent_value)

    print("The point of tangency on the stability boundary is:")
    print(f"z_t = {re_z:.4f} + {im_z:.4f}i")
    print("\nThe stability angle alpha is given by the equation:")
    print(f"alpha = arctan(| Im(z_t) / Re(z_t) |)")
    print(f"alpha = arctan(| {im_z:.4f} / {re_z:.4f} |)")
    print(f"alpha = arctan({tangent_value:.4f})")
    
    # For the exact expression, tangent_value is 3 * sqrt(3)
    exact_tangent_expr = "3 * sqrt(3)"
    print(f"\nThe exact expression is alpha = arctan({exact_tangent_expr})")

    print(f"\nThe numerical value of alpha is {alpha_rad:.4f} radians.")

solve_bdf4_angle()