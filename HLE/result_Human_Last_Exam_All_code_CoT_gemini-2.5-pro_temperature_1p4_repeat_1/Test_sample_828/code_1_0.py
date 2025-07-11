import sympy

def solve_tensor_transformation():
    """
    This function performs a symbolic transformation of the given stress-energy tensor
    from Cartesian to polar coordinates to find the component T_thetatheta and solve for K.
    """

    # 1. Define symbolic variables from the problem
    T_cal, a, omega, r, theta = sympy.symbols('mathcal{T} a omega r theta')

    # 2. Define the transformation derivatives from x=r*cos(theta), y=r*sin(theta)
    dx_dtheta = -r * sympy.sin(theta)
    dy_dtheta = r * sympy.cos(theta)

    # 3. Define the given Stress-Energy tensor components in Cartesian coordinates
    # We use a helper variable C for the common factor to keep expressions tidy.
    C = T_cal * (a * omega)**2
    T_xx = C * sympy.sin(theta)**2
    T_yy = C * sympy.cos(theta)**2
    T_xy = -C * sympy.sin(theta) * sympy.cos(theta)

    # 4. Apply the transformation law to find T'_thetatheta
    # T'_thetatheta = (dx/dtheta)^2*T_xx + (dy/dtheta)^2*T_yy + 2*(dx/dtheta)*(dy/dtheta)*T_xy
    T_thetatheta_polar = (dx_dtheta)**2 * T_xx + \
                       (dy_dtheta)**2 * T_yy + \
                       2 * dx_dtheta * dy_dtheta * T_xy

    # 5. Evaluate on the ring of radius 'a' by substituting r=a
    T_thetatheta_polar_at_a = T_thetatheta_polar.subs(r, a)

    # 6. Simplify the expression for T'_thetatheta
    T_thetatheta_simplified = sympy.simplify(T_thetatheta_polar_at_a)
    
    print("The transformed tensor component T'_thetatheta in polar coordinates is:")
    print(f"T'_thetatheta = {T_thetatheta_simplified}\n")

    # 7. Solve for K using the relation: T'_thetatheta = a^2*sin^2(theta)*T + K
    # K = T'_thetatheta - a^2*sin^2(theta)*T
    K_expr = T_thetatheta_simplified - a**2 * sympy.sin(theta)**2 * T_cal
    
    # Factor the expression for K for a more compact representation
    K_factored = sympy.factor(K_expr)

    # 8. Print the final equation for K
    # The output represents the equation: K = T*a^2*(a^2*w^2 - sin^2(theta))
    # where T is mathcal{T} and w is omega.
    print("From the relation T'_thetatheta = a^2*sin^2(theta)*T + K, we solve for K:")
    print("K = T'_thetatheta - a^2*sin^2(theta)*T")
    print("\nSubstituting the calculated value of T'_thetatheta gives the final expression for K:")
    print(f"K = {K_factored}")

if __name__ == '__main__':
    solve_tensor_transformation()