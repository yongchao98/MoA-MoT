def derive_critical_speed():
    """
    This function programmatically derives and displays the formula for the critical
    speed of an oversteering vehicle based on the linear single-track model.
    """
    
    # --- Step 1: Define parameters and state the initial equation ---
    print("Derivation of Critical Speed (v_crit)")
    print("------------------------------------------")
    print("This derivation uses the following parameters of the linear single-track model:")
    print("  a: distance from Center of Gravity (CG) to front axle")
    print("  b: distance from CG to rear axle")
    print("  c_f: cornering stiffness of the front axle")
    print("  c_r: cornering stiffness of the rear axle")
    print("  m: vehicle mass")
    print("  I: vehicle moment of inertia about the yaw axis")
    print("  v: constant forward speed")
    print("\n")

    # --- Step 2: Explain the stability condition ---
    print("Step 1: The Stability Condition")
    print("The stability of the vehicle's lateral dynamics is determined by the sign of the constant term 'Q' in the system's characteristic equation.")
    print("The system is stable if Q > 0. The onset of instability (critical speed) occurs when Q = 0.")
    print("The expression for Q is:")
    print("Q = ( (a + b)^2 * c_f * c_r ) / ( m * I * v^2 ) - ( a * c_f - b * c_r ) / I")
    print("\n")

    # --- Step 3: Set up the equation for critical speed ---
    print("Step 2: Setting up the equation for v_crit")
    print("We set Q = 0 to find the critical speed v_crit:")
    print("( (a + b)^2 * c_f * c_r ) / ( m * I * v_crit^2 ) - ( a * c_f - b * c_r ) / I = 0")
    print("\n")
    
    # --- Step 4: Rearrange the equation ---
    print("Step 3: Rearranging the equation to solve for v_crit^2")
    print("We move one term to the other side of the equation:")
    print("( (a + b)^2 * c_f * c_r ) / ( m * I * v_crit^2 ) = ( a * c_f - b * c_r ) / I")
    print("Note: For an oversteering vehicle, the term (a*c_f - b*c_r) is positive, ensuring a real solution exists.")
    print("\n")

    # --- Step 5: Solve for v_crit^2 ---
    print("Step 4: Isolating v_crit^2")
    print("By cross-multiplication and simplification (canceling 'I' from both sides), we get:")
    print("v_crit^2 = ( (a + b)^2 * c_f * c_r ) / ( m * (a * c_f - b * c_r) )")
    print("\n")

    # --- Step 6: Final Formula ---
    print("Step 5: Final Formula for Critical Speed")
    print("Taking the square root of both sides gives the final formula for the critical speed:")
    print("v_crit = sqrt( [ (a + b)^2 * c_f * c_r ] / [ m * (a * c_f - b * c_r) ] )")
    print("\n")

    # --- Displaying the final equation with each component ---
    print("The final equation is composed of the following parts:")
    print("Numerator = (a + b)^2 * c_f * c_r")
    print("Denominator = m * (a * c_f - b * c_r)")
    print("v_crit = sqrt( Numerator / Denominator )")


if __name__ == '__main__':
    derive_critical_speed()
