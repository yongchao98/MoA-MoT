import sympy
from sympy import symbols, Matrix, sin, cos, diff, Eq, solve, factor

def solve_tensor_transformation():
    """
    This function performs a coordinate transformation on a stress-energy tensor
    and solves for the unknown component K.
    """
    # 1. Define all symbolic variables
    # Tcal is used for \mathcal{T} to avoid conflict with SymPy's matrix transpose .T
    Tcal, a, omega, r, theta = symbols('Tcal a omega r theta')

    # 2. Define the transformation from polar to Cartesian coordinates for x and y
    # x^1 = x, x^2 = y. x'^2 = r, x'^3 = theta
    # In the problem, we evaluate the tensor on the ring, so r = a.
    x_coord = a * cos(theta)
    y_coord = a * sin(theta)

    # 3. Calculate the necessary partial derivatives for the transformation matrix
    # We need the derivatives of Cartesian coordinates (x,y) w.r.t. the polar coordinate theta.
    # Note: ∂x/∂θ = -a*sin(θ), ∂y/∂θ = a*cos(θ)
    # The problem implies we are at a fixed radius a, so r is replaced by a.
    # We are calculating T'_θθ, so we only need derivatives w.r.t theta.
    # ∂x^μ/∂x'^θ for μ=x,y. Others are 0.
    dx_dtheta = diff(a * cos(theta), theta)
    dy_dtheta = diff(a * sin(theta), theta)

    # 4. Define the non-zero Cartesian stress-energy tensor components (spatial part)
    # The given tensor components depend on the polar angle theta, which is unusual but we follow the prompt.
    S = (a * omega)**2
    T_xx = Tcal * S * sin(theta)**2
    T_yy = Tcal * S * cos(theta)**2
    T_xy = -Tcal * S * sin(theta) * cos(theta)
    T_yx = T_xy # The tensor is symmetric

    # 5. Apply the tensor transformation law to find T'_θθ
    # T'_θθ = (∂x/∂θ)² T_xx + (∂y/∂θ)² T_yy + 2(∂x/∂θ)(∂y/∂θ) T_xy
    T_thetatheta = (dx_dtheta**2 * T_xx +
                    dy_dtheta**2 * T_yy +
                    2 * dx_dtheta * dy_dtheta * T_xy)

    # 6. Simplify the resulting expression
    T_thetatheta_simplified = sympy.simplify(T_thetatheta)
    
    # 7. Solve for K using the equation provided in the problem
    # T'_θθ = a² * sin²θ * Tcal + K
    K = symbols('K')
    equation = Eq(T_thetatheta_simplified, a**2 * sin(theta)**2 * Tcal + K)
    
    # Solve for K
    solution = solve(equation, K)
    K_expression = solution[0]
    
    # Factor the expression for a cleaner look
    K_factored = factor(K_expression)
    
    # 8. Print the results. The final equation is K = expression.
    # We rename Tcal back to T for the final output to match the problem statement.
    final_K = str(K_factored).replace('Tcal', 'T')
    
    print("The transformed tensor component T_thetatheta is:")
    print(str(T_thetatheta_simplified).replace('Tcal', 'T'))
    print("\nThe equation given is T_thetatheta = T*a**2*sin(theta)**2 + K")
    print("Solving for K, we get:")
    print(f"K = {final_K}")

solve_tensor_transformation()