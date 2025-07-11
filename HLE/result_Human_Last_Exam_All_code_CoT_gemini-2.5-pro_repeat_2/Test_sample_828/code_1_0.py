import sympy as sp

def solve_for_K():
    """
    This function calculates the expression for K by first transforming the
    stress-energy tensor to polar coordinates and then solving the given equation.
    """
    # 1. Define symbolic variables
    # T_cal represents the script T in the problem
    T_cal, a, omega, theta = sp.symbols('T_cal a omega theta')

    # 2. Define the Cartesian components of the stress-energy tensor T_munu
    # We only need the spatial components for this transformation.
    T_xx = T_cal * (a * omega)**2 * sp.sin(theta)**2
    T_xy = -T_cal * (a * omega)**2 * sp.sin(theta) * sp.cos(theta)
    T_yx = T_xy
    T_yy = T_cal * (a * omega)**2 * sp.cos(theta)**2

    # 3. Define the partial derivatives for the coordinate transformation x=r*cos(theta), y=r*sin(theta)
    # The problem is defined on a ring of radius 'a', so we set r=a.
    dx_dtheta = -a * sp.sin(theta)
    dy_dtheta = a * sp.cos(theta)

    # 4. Calculate the T_thetatheta component using the tensor transformation law
    # T'_ab = (dx^c/dx'^a) * (dx^d/dx'^b) * T_cd
    T_thetatheta = (dx_dtheta**2 * T_xx) + \
                  (dx_dtheta * dy_dtheta * T_xy) + \
                  (dy_dtheta * dx_dtheta * T_yx) + \
                  (dy_dtheta**2 * T_yy)

    # Simplify the expression for T_thetatheta
    T_thetatheta_simplified = sp.simplify(T_thetatheta)

    print(f"The calculated expression for T_thetatheta is: {T_thetatheta_simplified}")
    
    # 5. Define the given term from the problem statement
    given_term = a**2 * sp.sin(theta)**2 * T_cal

    # 6. Solve for K using the equation: T_thetatheta = given_term + K
    K = T_thetatheta_simplified - given_term
    
    # Factor the expression for K to make it cleaner
    K_simplified = sp.factor(K)
    
    print("\nThe problem states: T_thetatheta = a**2 * sin(theta)**2 * T_cal + K")
    print(f"Solving for K, we get:")
    print(f"K = {K_simplified}")

# Execute the function
solve_for_K()