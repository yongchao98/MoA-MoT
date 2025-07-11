import sympy

def solve_tensor_transformation():
    """
    This function calculates the T_thetatheta component of the stress-energy tensor
    in polar coordinates and then determines the value of K based on the expression
    given in the problem.
    """
    # Define the symbols used in the problem
    T, a, w, r, theta = sympy.symbols('mathcal{T} a omega r theta')

    # Cartesian components of the stress-energy tensor are given.
    # We interpret the 'theta' in these expressions as the polar coordinate angle.
    T_xx = T * (a*w)**2 * sympy.sin(theta)**2
    T_xy = -T * (a*w)**2 * sympy.sin(theta) * sympy.cos(theta)
    T_yy = T * (a*w)**2 * sympy.cos(theta)**2

    # Partial derivatives for the coordinate transformation x = r*cos(theta), y = r*sin(theta)
    # with respect to theta.
    dx_dtheta = sympy.diff(r * sympy.cos(theta), theta)
    dy_dtheta = sympy.diff(r * sympy.sin(theta), theta)

    # Use the tensor transformation law to find the T_thetatheta component in polar coordinates.
    # T'_{thetatheta} = (dx/dtheta)^2 * T_xx + 2*(dx/dtheta)*(dy/dtheta)*T_xy + (dy/dtheta)^2 * T_yy
    T_thetatheta_transformed = sympy.simplify(
        dx_dtheta**2 * T_xx +
        2 * dx_dtheta * dy_dtheta * T_xy +
        dy_dtheta**2 * T_yy
    )

    # The problem describes a ring of radius 'a', so we evaluate the component at r = a.
    T_thetatheta_at_a = T_thetatheta_transformed.subs(r, a)

    # The problem states that T_thetatheta is of the form: a^2 * sin^2(theta) * T + K
    # We set our calculated result equal to this form and solve for K.
    # K = T_thetatheta_at_a - a^2*sin^2(theta)*T
    K = sympy.simplify(T_thetatheta_at_a - a**2 * sympy.sin(theta)**2 * T)

    # The problem asks to output the final equation. Let's present the derived components clearly.
    # The derived T_thetatheta is T_thetatheta_at_a.
    # The equation form is T_thetatheta = a**2 * sympy.sin(theta)**2 * T + K
    print(f"The transformed component is: T_thetatheta = {T_thetatheta_at_a}")
    print(f"Given the form: T_thetatheta = a**2*sin(theta)**2*T + K")
    print(f"Therefore, K is: {K}")
    
    # We will now construct the final equation as requested by the prompt.
    term1 = sympy.pretty(a**2 * sympy.sin(theta)**2 * T, use_unicode=True)
    k_term = sympy.pretty(K, use_unicode=True)
    
    # Output each number in the final equation.
    print("\nThe final equation with K substituted is:")
    print(f"T_θθ = {term1} + ({k_term})")
    

solve_tensor_transformation()

# The question is 'what is K?'. We extract the expression for K as the final answer.
T, a, w, theta = sympy.symbols('mathcal{T} a omega theta')
K = T*a**2*((a*w)**2 - sympy.sin(theta)**2)
print(f"\n<<<{sympy.sstr(K, full_obj=False)}>>>")
