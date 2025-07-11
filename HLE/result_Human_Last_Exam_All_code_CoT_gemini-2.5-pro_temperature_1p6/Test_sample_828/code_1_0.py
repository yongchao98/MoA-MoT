import sympy

def solve_tensor_problem():
    """
    This function calculates the T_thetatheta component of the stress-energy tensor
    in polar coordinates and then determines the value of K based on the problem statement.
    """

    # Step 1: Define all symbolic variables
    a, omega, r, T_cal = sympy.symbols('a omega r T_cal')
    theta = sympy.Symbol('theta')

    # Step 2: Define the transformation from polar to Cartesian coordinates
    # The problem asks to transform from Cartesian to Polar, so we use the inverse relations
    # x = r*cos(theta), y = r*sin(theta)
    
    # Step 3: Define the non-zero spatial components of the Cartesian stress-energy tensor T_mu_nu
    # The factor T_cal corresponds to script T in the problem
    T_xx = T_cal * (a * omega)**2 * sympy.sin(theta)**2
    T_yy = T_cal * (a * omega)**2 * sympy.cos(theta)**2
    T_xy = -T_cal * (a * omega)**2 * sympy.sin(theta) * sympy.cos(theta)
    T_yx = T_xy

    # Step 4: Calculate the necessary partial derivatives for the transformation
    # We need d(x)/d(theta) and d(y)/d(theta)
    x = r * sympy.cos(theta)
    y = r * sympy.sin(theta)
    
    dx_dtheta = sympy.diff(x, theta)
    dy_dtheta = sympy.diff(y, theta)

    # Step 5: The transformation law for the T_thetatheta component is:
    # T'_thetatheta = (dx/dtheta)^2 * T_xx + (dy/dtheta)^2 * T_yy + 2 * (dx/dtheta)*(dy/dtheta) * T_xy
    T_thetatheta_prime = (dx_dtheta**2 * T_xx + 
                          dy_dtheta**2 * T_yy + 
                          2 * dx_dtheta * dy_dtheta * T_xy)

    # The dust ring has radius 'a', so we substitute r=a
    T_thetatheta_prime_on_ring = T_thetatheta_prime.subs(r, a)

    # Step 6: Simplify the expression
    T_thetatheta_simplified = sympy.simplify(T_thetatheta_prime_on_ring)

    # Print out the derivation steps
    print("Step 1: The original Cartesian tensor components are:")
    print(f"T_xx = {T_xx}")
    print(f"T_yy = {T_yy}")
    print(f"T_xy = {T_xy}")
    print("-" * 20)
    print("Step 2: The partial derivatives for the transformation are:")
    print(f"dx/dtheta = {dx_dtheta}")
    print(f"dy/dtheta = {dy_dtheta}")
    print("-" * 20)
    print("Step 3: The transformation formula for T'_thetatheta at r=a is:")
    # Print the equation with substituted parts. Using a temp variable for printing clarity
    term1 = f"({dx_dtheta.subs(r, a)})**2 * ({T_xx})"
    term2 = f"({dy_dtheta.subs(r, a)})**2 * ({T_yy})"
    term3 = f"2 * ({dx_dtheta.subs(r, a)}) * ({dy_dtheta.subs(r, a)}) * ({T_xy})"
    print(f"T'_thetatheta = {term1} + {term2} + {term3}")
    print("-" * 20)
    print("Step 4: The simplified result for T'_thetatheta is:")
    print(f"T'_thetatheta = {T_thetatheta_simplified}")
    print("-" * 20)

    # Step 7: The problem states T'_thetatheta = a^2*sin(theta)^2*T_cal + K
    # Equating our result with the given expression:
    K = sympy.Symbol('K')
    given_form = a**2 * T_cal * sympy.sin(theta)**2 + K
    
    equation_for_K = sympy.Eq(T_thetatheta_simplified, given_form)
    
    print("Step 5: Equating the derived result with the form given in the problem:")
    print(f"{T_thetatheta_simplified} = {given_form}")
    print("-" * 20)

    # Step 8: Solve for K and analyze the result
    K_solution = sympy.solve(equation_for_K, K)[0]
    print(f"Step 6: Solving for K gives:")
    print(f"K = {K_solution}")
    print("-" * 20)
    
    print("Step 7: For K to be a constant, it cannot depend on theta.")
    print("The expression for K depends on theta unless the coefficient of the sin(theta)**2 term is zero.")
    print("The term is -a**2*T_cal*sin(theta)**2. For this to be zero for all theta, T_cal must be 0.")
    print("If T_cal = 0, the original tensor is all zeros.")
    print("Then T'_thetatheta is also 0.")
    final_K = K_solution.subs(T_cal, 0)
    print("Substituting T_cal = 0 into the expression for K, we get:")
    print(f"K = {final_K}")

    return final_K

if __name__ == "__main__":
    final_answer = solve_tensor_problem()
    # The final answer is requested in a specific format
    # print(f"<<<{final_answer}>>>")
    
solve_tensor_problem()
<<<0>>>