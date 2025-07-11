def derive_critical_speed():
    """
    This function prints the step-by-step derivation for the critical speed
    of an oversteering vehicle based on the linear single-track model.
    """
    print("### Derivation of Critical Speed for an Oversteering Vehicle ###\n")

    print("Step 1: State-Space Representation of the Single-Track Model")
    print("------------------------------------------------------------")
    print("The vehicle's linear lateral dynamics are described by a state-space model d(x)/dt = A * x,")
    print("where the state vector x = [β, r]^T consists of the sideslip angle (β) and yaw rate (r).")
    print("The state matrix A is given by:")
    print("A = [[ -(c_f + c_r)/(m*v),  -(a*c_f - b*c_r)/(m*v^2) - 1 ],")
    print("     [ -(a*c_f - b*c_r)/I, -(a^2*c_f + b^2*c_r)/(I*v)    ]]\n")

    print("Step 2: Condition for Dynamic Stability")
    print("--------------------------------------")
    print("For the vehicle to be stable, the real parts of all eigenvalues of matrix A must be negative.")
    print("This requires two conditions to be met:")
    print("  a) The trace of the matrix, trace(A), must be negative.")
    print("  b) The determinant of the matrix, det(A), must be positive.")
    print("The trace is always negative for a forward speed v > 0, so stability depends on the determinant.\n")

    print("Step 3: Calculating the Determinant")
    print("-----------------------------------")
    print("The determinant of A is det(A) = A[0,0]*A[1,1] - A[0,1]*A[1,0].")
    print("After algebraic simplification, the expression for the determinant is:")
    print("det(A) = (c_f * c_r * (a+b)^2) / (m * I * v^2) - (a*c_f - b*c_r) / I\n")

    print("Step 4: The Oversteering Condition and Critical Speed")
    print("------------------------------------------------------")
    print("A vehicle is defined as 'oversteering' if the term (a*c_f - b*c_r) is positive.")
    print("For an oversteering vehicle, the stability condition det(A) > 0 can be violated as speed 'v' increases.")
    print("The 'critical speed' (v_crit) is the speed at which stability is lost, which occurs when det(A) = 0.\n")

    print("Step 5: Solving for the Critical Speed (v_crit)")
    print("------------------------------------------------")
    print("We set the determinant to zero and solve for v = v_crit:")
    print("\nEquation to solve:")
    print("(c_f * c_r * (a + b)^2) / (m * I * v_crit^2) - (a * c_f - b * c_r) / I = 0")
    
    print("\nRearranging the equation:")
    print("(c_f * c_r * (a + b)^2) / (m * I * v_crit^2) = (a * c_f - b * c_r) / I")
    
    print("\nSolving for v_crit^2:")
    print("v_crit^2 = (c_f * c_r * (a + b)^2 * I) / (m * I * (a * c_f - b * c_r))")
    print("v_crit^2 = (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r))\n")

    print("### Final Result ###")
    print("Taking the square root gives the final expression for the critical speed:")
    print("\nv_crit = sqrt( (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r)) )")


if __name__ == '__main__':
    derive_critical_speed()
    final_formula = "sqrt( (c_f * c_r * (a + b)**2) / (m * (a * c_f - b * c_r)) )"
    # The final answer is the formula itself.
    # No numerical value can be computed as inputs are not provided.
    # Let's consider the formula as the answer content.
    print(f"\n<<<v_crit = {final_formula}>>>")