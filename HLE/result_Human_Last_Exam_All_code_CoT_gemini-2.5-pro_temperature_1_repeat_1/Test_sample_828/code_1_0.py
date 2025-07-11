import sympy as sp

def solve_tensor_transformation():
    """
    This function calculates the transformed stress-energy tensor component T'_{\u03B8\u03B8}
    and then solves for K based on the equation provided in the problem.
    """
    # Define the symbolic variables
    # r: polar radius
    # theta: polar angle
    # a: radius of the dust ring
    # omega: angular velocity
    # calT: a constant related to dust density and mass (mathcal{T})
    r, theta, a, omega = sp.symbols('r theta a omega', real=True, positive=True)
    calT = sp.Symbol('calT')

    # Define the Cartesian coordinates in terms of polar coordinates
    x = r * sp.cos(theta)
    y = r * sp.sin(theta)

    # Given components of the stress-energy tensor in Cartesian coordinates
    # Let V2 = (a*omega)^2 for simplicity
    V2 = (a * omega)**2
    T_xx = calT * V2 * sp.sin(theta)**2
    T_yy = calT * V2 * sp.cos(theta)**2
    T_xy = -calT * V2 * sp.sin(theta) * sp.cos(theta)

    # Calculate the necessary partial derivatives for the transformation
    # We need dx/d\theta and dy/d\theta
    dx_dtheta = sp.diff(x, theta)
    dy_dtheta = sp.diff(y, theta)

    # Apply the tensor transformation law for the T'_{\u03B8\u03B8} component
    # T'_{\u03B8\u03B8} = (dx/d\theta)^2 * T_xx + (dy/d\theta)^2 * T_yy + 2*(dx/d\theta)*(dy/d\theta)*T_xy
    T_prime_thetatheta = (dx_dtheta**2 * T_xx +
                          dy_dtheta**2 * T_yy +
                          2 * dx_dtheta * dy_dtheta * T_xy)

    # Simplify the resulting expression
    T_prime_thetatheta_simplified = sp.simplify(T_prime_thetatheta)

    # The problem specifies a ring of radius 'a', so we evaluate our result at r=a
    T_prime_thetatheta_at_a = T_prime_thetatheta_simplified.subs(r, a)

    # The problem states: T'_{\u03B8\u03B8} = a^2 * sin^2(\theta) * calT + K
    # We now solve for K
    given_term = a**2 * sp.sin(theta)**2 * calT
    K = T_prime_thetatheta_at_a - given_term
    K_simplified = sp.simplify(K)

    # Print the final equation with each term, as requested.
    # We use unicode characters for a clear mathematical representation.
    print("Step 1: The calculated transformed component T'(\u03B8,\u03B8) at r=a is:")
    print(f"T'(\u03B8,\u03B8) = {sp.pretty(T_prime_thetatheta_at_a, use_unicode=True)}")
    print("\nStep 2: The problem provides the equation relating this component to K:")
    print(f"{sp.pretty(T_prime_thetatheta_at_a, use_unicode=True)} = {sp.pretty(given_term, use_unicode=True)} + K")
    print("\nStep 3: Solving for K, we get:")
    print(f"K = {sp.pretty(T_prime_thetatheta_at_a, use_unicode=True)} - ({sp.pretty(given_term, use_unicode=True)})")
    print("\nStep 4: The final simplified expression for K is:")
    print(f"K = {sp.pretty(K_simplified, use_unicode=True)}")

if __name__ == '__main__':
    solve_tensor_transformation()