import sympy

def solve_tensor_transformation():
    """
    This function performs the symbolic transformation of the stress-energy tensor
    from Cartesian to polar coordinates to find the component T_{\theta\theta} and then solves for K.
    """
    # 1. Define symbolic variables
    # Using 'T_cal' for the script T
    T_cal, a, omega, theta = sympy.symbols('T_cal a omega theta')

    # 2. Define the Cartesian tensor components (on the ring r=a)
    # The problem gives T_{mu,nu} = T_cal * Matrix
    # So we extract the spatial components we need
    T_xx = T_cal * (a * omega)**2 * sympy.sin(theta)**2
    T_yy = T_cal * (a * omega)**2 * sympy.cos(theta)**2
    T_xy = -T_cal * (a * omega)**2 * sympy.sin(theta) * sympy.cos(theta)
    # T_yx is the same as T_xy since the tensor is symmetric
    T_yx = T_xy

    # 3. Calculate the partial derivatives for the transformation x=r*cos(theta), y=r*sin(theta)
    # We evaluate at the ring's radius r=a
    # dx/d(theta) evaluated at r=a
    dx_dtheta = -a * sympy.sin(theta)
    # dy/d(theta) evaluated at r=a
    dy_dtheta = a * sympy.cos(theta)

    # 4. Apply the tensor transformation rule for T'_{\theta\theta}
    # T'_{\theta\theta} = (dx/d\theta)^2 * T_xx + (dy/d\theta)^2 * T_yy + 2*(dx/d\theta)*(dy/d\theta)*T_xy
    T_prime_thetatheta = (dx_dtheta**2 * T_xx) + \
                          (dy_dtheta**2 * T_yy) + \
                          (dx_dtheta * dy_dtheta * T_xy) + \
                          (dy_dtheta * dx_dtheta * T_yx)

    # 5. Simplify the expression for T'_{\theta\theta}
    T_prime_thetatheta_simplified = sympy.simplify(T_prime_thetatheta)

    # 6. Set up the equation from the problem statement to solve for K
    K = sympy.Symbol('K')
    # The equation is T'_{\theta\theta} = a^2 * sin^2(theta) * T_cal + K
    equation_for_K = sympy.Eq(T_prime_thetatheta_simplified, a**2 * sympy.sin(theta)**2 * T_cal + K)

    # 7. Solve for K
    K_solution = sympy.solve(equation_for_K, K)

    # The result K_solution is a list containing one element, which is our expression for K
    K_expression = K_solution[0]
    
    # We can factor the expression for a cleaner look
    K_expression_factored = sympy.factor(K_expression)

    # 8. Print the result in the format of an equation.
    # We will output the factored form.
    final_equation = sympy.Eq(K, K_expression_factored)
    print("The transformation yields the component in polar coordinates:")
    print(sympy.Eq(sympy.Symbol("T_thetatheta"), T_prime_thetatheta_simplified))
    print("\nGiven the relation from the problem:")
    print(f"T_thetatheta = a^2*sin(theta)^2*T_cal + K")
    print("\nSolving for K, we get:")
    # We construct and print the final equation for K
    print(f"{final_equation.lhs} = {final_equation.rhs}")


if __name__ == "__main__":
    solve_tensor_transformation()