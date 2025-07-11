import sympy

def derive_critical_speed():
    """
    Derives the critical speed for an oversteering vehicle using the linear single-track model.
    """
    # 1. Define symbolic variables for the vehicle parameters
    a, b, c_f, c_r, m, I, v = sympy.symbols('a b c_f c_r m I v', positive=True)

    print("Step 1: Define the state-space matrix A for the linear single-track model.")
    print("The state vector is x = [beta, r]^T, where beta is the sideslip angle and r is the yaw rate.")
    
    # 2. Define the elements of the state matrix A
    # A = [[A11, A12], [A21, A22]]
    # These elements are derived from the vehicle's equations of motion:
    # m*v*(d(beta)/dt + r) = -(c_f + c_r)*beta - (a*c_f - b*c_r)/v * r
    # I*d(r)/dt = -(a*c_f - b*c_r)*beta - (a**2*c_f + b**2*c_r)/v * r
    
    A11 = -(c_f + c_r) / (m * v)
    A12 = (b * c_r - a * c_f) / (m * v**2) - 1
    A21 = (b * c_r - a * c_f) / I
    A22 = -(a**2 * c_f + b**2 * c_r) / (I * v)

    A = sympy.Matrix([[A11, A12], [A21, A22]])
    
    print("\nThe state matrix A is:")
    sympy.pprint(A)

    print("\nStep 2: Calculate the determinant of the state matrix A.")
    # 3. Calculate the determinant of A
    det_A = A.det()
    # The raw determinant is complex, so we simplify it.
    det_A_simplified = sympy.simplify(det_A)
    
    print("\nThe determinant is:")
    sympy.pprint(det_A_simplified)
    
    print("\nNote: The term (a*c_f - b*c_r) defines the vehicle's handling characteristics.")
    print("For an OVERSTEERING vehicle, a*c_f > b*c_r.")
    print("For an UNDERSTEERING vehicle, a*c_f < b*c_r.")
    print("\nStability requires det(A) > 0. For an oversteering vehicle, the term (b*c_r - a*c_f) is negative.")
    print("As speed 'v' increases, the first term in the determinant decreases, and det(A) can become zero or negative, leading to instability.")

    print("\nStep 3: Find the critical speed by setting the determinant to zero and solving for v.")
    
    # 4. Set the determinant to zero to find the critical speed
    # We solve for v**2 first
    v_crit_squared_eq = sympy.Eq(det_A_simplified, 0)
    
    # The solve function returns a list of solutions
    v_crit_squared_sol = sympy.solve(v_crit_squared_eq, v**2)
    
    # We expect one positive solution for v**2
    v_crit_squared = v_crit_squared_sol[0]
    
    # Let's clean up the expression for better readability
    # Numerator: (a+b)^2 * c_f * c_r
    # Denominator: m * (a*c_f - b*c_r)
    num = (a + b)**2 * c_f * c_r
    den = m * (a * c_f - b * c_r)
    v_crit_squared_final = num / den

    print("\nThe equation for the square of the critical speed (v_crit^2) is:")
    print("v_crit^2 = ", end="")
    sympy.pprint(v_crit_squared_final)

    # 5. Take the square root to find the final expression for v_crit
    v_crit = sympy.sqrt(v_crit_squared_final)

    print("\nStep 4: The final derived formula for the critical speed (v_crit) is:")
    print("\nv_crit = sqrt( (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r)) )")
    print("\nSymbolic representation:")
    print("v_crit = ", end="")
    sympy.pprint(v_crit)
    
    print("\nNote that the vehicle's moment of inertia 'I' cancels out and does not appear in the final formula for critical speed.")

if __name__ == '__main__':
    derive_critical_speed()
