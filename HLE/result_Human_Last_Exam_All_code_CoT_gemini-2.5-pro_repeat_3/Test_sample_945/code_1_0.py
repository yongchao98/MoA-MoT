def derive_critical_speed():
    """
    This function prints the step-by-step derivation of the critical speed
    for an oversteering vehicle using the linear single-track model.
    """
    
    # Plan of the derivation:
    # 1. State the parameters and the equations of motion in state-space form.
    # 2. Define the state matrix A for the system.
    # 3. Explain the stability condition based on the characteristic equation of the system.
    #    Instability occurs when the constant term of the polynomial, C2 = det(A), becomes non-positive.
    # 4. Calculate the determinant of the state matrix A.
    # 5. Set the determinant to zero to find the critical speed v_crit where stability is lost.
    # 6. Present the final formula for v_crit.
    
    print("Derivation of Critical Speed for an Oversteering Vehicle\n")
    print("="*60)

    # Step 1: Define parameters and equations of motion
    print("Step 1: Parameters and Equations of Motion")
    print("-" * 40)
    print("The linear single-track model uses the following parameters:")
    print("  a:   distance from Center of Gravity (CG) to the front axle")
    print("  b:   distance from CG to the rear axle")
    print("  c_f: cornering stiffness of the front axle")
    print("  c_r: cornering stiffness of the rear axle")
    print("  m:   vehicle mass")
    print("  I:   vehicle moment of inertia about the vertical axis")
    print("  v:   vehicle longitudinal speed (assumed constant)\n")
    
    print("The system's state is described by the side-slip angle (beta) and yaw rate (r).")
    print("The linearized equations of motion (x_dot = A*x) for a steering angle of zero are:")
    print("  beta_dot = -((c_f + c_r)/(m*v)) * beta - (1 + (a*c_f - b*c_r)/(m*v^2)) * r")
    print("  r_dot    = -((a*c_f - b*c_r)/I) * beta - ((a^2*c_f + b^2*c_r)/(I*v)) * r\n")

    # Step 2: Define the State Matrix A
    print("Step 2: State-Space Matrix A")
    print("-" * 40)
    print("From the equations above, the state matrix A is:")
    print("  A = [ -((c_f + c_r)/(m*v))      -(1 + (a*c_f - b*c_r)/(m*v^2)) ]")
    print("      [ -((a*c_f - b*c_r)/I)       -((a^2*c_f + b^2*c_r)/(I*v))     ]\n")

    # Step 3: Stability Condition
    print("Step 3: Stability Condition")
    print("-" * 40)
    print("The system's stability depends on the eigenvalues of matrix A.")
    print("The eigenvalues are the roots of the characteristic equation: lambda^2 + C1*lambda + C2 = 0.")
    print("For stability, all coefficients (C1, C2) must be positive.")
    print("The system becomes unstable when C2 <= 0. The critical speed is the speed 'v' at which C2 = 0.")
    print("The coefficient C2 is the determinant of matrix A: C2 = det(A).\n")

    # Step 4: Calculate det(A)
    print("Step 4: Calculating the Determinant of A")
    print("-" * 40)
    print("det(A) = [(-c_f - c_r)/(m*v)] * [(-a^2*c_f - b^2*c_r)/(I*v)] - [-(1 + (a*c_f - b*c_r)/(m*v^2))] * [-(a*c_f - b*c_r)/I]")
    print("\nSimplifying the expression gives:")
    print("det(A) = [(c_f + c_r)*(a^2*c_f + b^2*c_r)] / [m*I*v^2] - [(a*c_f - b*c_r)/I] * [1 + (a*c_f - b*c_r)/(m*v^2)]")
    print("\nFurther simplification by expanding and collecting terms over a common denominator leads to:")
    print("det(A) = [c_f*c_r*(a+b)^2 - m*v^2*(a*c_f - b*c_r)] / [m*I*v^2]\n")
    
    # Step 5: Solve for Critical Speed
    print("Step 5: Solving for Critical Speed (v_crit)")
    print("-" * 40)
    print("To find the critical speed, we set det(A) = 0:")
    print("  [c_f*c_r*(a+b)^2 - m*v_crit^2*(a*c_f - b*c_r)] / [m*I*v_crit^2] = 0")
    print("\nThis requires the numerator to be zero:")
    print("  c_f*c_r*(a+b)^2 - m*v_crit^2*(a*c_f - b*c_r) = 0")
    print("\nRearranging to solve for v_crit^2:")
    print("  m*v_crit^2*(a*c_f - b*c_r) = c_f*c_r*(a+b)^2")
    print("  v_crit^2 = (c_f * c_r * (a+b)^2) / (m * (a*c_f - b*c_r))\n")
    
    # Step 6: Final Formula
    print("Step 6: Final Formula for Critical Speed")
    print("-" * 40)
    print("An oversteering vehicle is defined by the condition (a*c_f - b*c_r) > 0.")
    print("This ensures that v_crit^2 is positive, yielding a real-valued critical speed.")
    print("Taking the square root gives the final formula:\n")
    
    print("Final Critical Speed Equation:")
    print("v_crit = sqrt( (c_f * c_r * (a + b)^2) / (m * (a*c_f - b*c_r)) )")
    print("="*60)

if __name__ == '__main__':
    derive_critical_speed()
