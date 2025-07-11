import sympy
from sympy import sin, cos, simplify, symbols, pretty_print

def solve_tensor_transformation():
    """
    This function calculates the T_thetatheta component of the stress-energy
    tensor in polar coordinates and then solves for K based on the problem statement.
    """
    # Define the symbolic variables used in the problem
    # T represents the calligraphic T in the problem
    a, omega, r, theta, T = symbols('a omega r theta T')

    # Cartesian coordinates in terms of polar coordinates
    x = r * cos(theta)
    y = r * sin(theta)

    # Calculate the required partial derivatives for the transformation
    dx_dtheta = sympy.diff(x, theta)
    dy_dtheta = sympy.diff(y, theta)

    # Define the given Cartesian stress-energy tensor components.
    # Let C = T * (a*omega)**2 for brevity in definition
    C = T * (a * omega)**2
    T_xx = C * sin(theta)**2
    T_yy = C * cos(theta)**2
    T_xy = -C * sin(theta) * cos(theta)

    # Apply the tensor transformation law to find T'_thetatheta
    # T'_thetatheta = (dx/dtheta)^2 * T_xx + (dy/dtheta)^2 * T_yy + 2*(dx/dtheta)*(dy/dtheta)*T_xy
    T_thetatheta_expr = dx_dtheta**2 * T_xx + dy_dtheta**2 * T_yy + 2 * dx_dtheta * dy_dtheta * T_xy

    # The problem considers a ring of radius 'a', so substitute r=a
    T_thetatheta_on_ring = T_thetatheta_expr.subs(r, a)

    # Simplify the resulting expression for T'_thetatheta
    T_thetatheta_simplified = simplify(T_thetatheta_on_ring)

    print("The calculated component T_thetatheta is:")
    pretty_print(T_thetatheta_simplified)
    print("-" * 20)
    
    # The problem asks to find K from the equation: T_thetatheta = a^2 * sin^2(theta) * T + K
    # We solve for K: K = T_thetatheta - a^2 * sin^2(theta) * T
    K_expr = T_thetatheta_simplified - T * a**2 * sin(theta)**2

    # Simplify the expression for K
    K_simplified = simplify(K_expr)

    # Print the final expression for K
    print("Given that T_thetatheta = a**2 * T * sin(theta)**2 + K,")
    print("The expression for K is:")
    print("K = ", end="")
    pretty_print(K_simplified)


if __name__ == '__main__':
    solve_tensor_transformation()
