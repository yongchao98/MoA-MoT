import sympy

def solve_for_K():
    """
    This function performs the tensor transformation and solves for K.
    """
    # Step 1: Define symbolic variables
    r, a, omega, theta = sympy.symbols('r a omega theta')
    T_cal = sympy.Symbol('T') # Representing the script T (mathcal{T})

    print("Step 1: Define symbolic variables for the problem.")
    print(f"Coordinates: r, theta")
    print(f"Constants: a, omega")
    print(f"Tensor density factor: {T_cal}")
    print("-" * 20)

    # Step 2: Define the Cartesian tensor components T_ij
    # Let C be the common factor T_cal * (a*omega)**2
    C = T_cal * (a * omega)**2
    T_xx = C * sympy.sin(theta)**2
    T_yy = C * sympy.cos(theta)**2
    T_xy = -C * sympy.sin(theta) * sympy.cos(theta)

    print("Step 2: Define the given Cartesian tensor components T_ij.")
    print(f"T_xx = {T_xx}")
    print(f"T_yy = {T_yy}")
    print(f"T_xy = {T_xy}")
    print("-" * 20)
    
    # Step 3: Define the partial derivatives for the coordinate transformation
    # x = r*cos(theta), y = r*sin(theta)
    dx_dtheta = -r * sympy.sin(theta)
    dy_dtheta = r * sympy.cos(theta)

    print("Step 3: Calculate the partial derivatives dx/dtheta and dy/dtheta.")
    print(f"dx/dtheta = {dx_dtheta}")
    print(f"dy/dtheta = {dy_dtheta}")
    print("-" * 20)
    
    # Step 4: Apply the tensor transformation law for T'_{\theta\theta}
    # T'_{\theta\theta} = (dx/dtheta)^2*T_xx + (dy/dtheta)^2*T_yy + 2*(dx/dtheta)*(dy/dtheta)*T_xy
    T_thetatheta_prime = (dx_dtheta**2 * T_xx +
                         dy_dtheta**2 * T_yy +
                         2 * dx_dtheta * dy_dtheta * T_xy)

    print("Step 4: Calculate the transformed component T'_{\\theta\\theta} using the transformation law.")
    # We print the un-simplified equation first to show the substitution
    print("Substituting into the formula T'_{\\theta\\theta} = (dx/dtheta)^2 * T_xx + (dy/dtheta)^2 * T_yy + 2*(dx/dtheta)*(dy/dtheta)*T_xy:")
    
    term1_str = f"({dx_dtheta})**2 * ({T_xx})"
    term2_str = f"({dy_dtheta})**2 * ({T_yy})"
    term3_str = f"2*({dx_dtheta})*({dy_dtheta}) * ({T_xy})"
    
    print(f"T'_{\\theta\\theta} = {term1_str} + {term2_str} + {term3_str}")
    
    # Simplify the expression
    T_thetatheta_prime_simplified = sympy.simplify(T_thetatheta_prime)
    
    print("\nSimplifying the expression for T'_{\\theta\\theta}:")
    print(f"T'_{\\theta\\theta} = {T_thetatheta_prime_simplified}")
    print("-" * 20)
    
    # Step 5: Evaluate at r=a for the ring
    T_thetatheta_final = T_thetatheta_prime_simplified.subs(r, a)
    
    print(f"Step 5: Evaluate T'_{\\theta\\theta} on the ring of radius r=a.")
    print(f"T_thetatheta(r=a) = {T_thetatheta_final}")
    print("-" * 20)

    # Step 6: Solve for K
    # Given: T_thetatheta = a**2 * sin(theta)**2 * T_cal + K
    # So, K = T_thetatheta_final - a**2 * T_cal * sympy.sin(theta)**2
    K = T_thetatheta_final - a**2 * T_cal * sympy.sin(theta)**2
    
    # We want to display the final result in the format K = ...
    K_simplified = sympy.simplify(K)
    
    print("Step 6: Solve for K from the equation T_thetatheta = a^2 * T * sin^2(theta) + K.")
    print(f"{T_thetatheta_final} = a**2*{T_cal}*sin(theta)**2 + K")
    print(f"K = {T_thetatheta_final} - a**2*{T_cal}*sin(theta)**2")
    
    print("\nFinal expression for K:")
    final_K_str = sympy.sstr(K_simplified, full_prec=True)
    # Make it more readable
    final_K_str = final_K_str.replace("omega", "w").replace("theta", "θ").replace("*", " * ")
    print(f"K = {final_K_str}")

# Run the solver
solve_for_K()