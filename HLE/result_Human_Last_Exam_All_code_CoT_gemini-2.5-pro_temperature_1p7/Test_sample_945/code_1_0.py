import math

def calculate_critical_speed():
    """
    Derives and calculates the critical speed for an oversteering vehicle.

    This function first prints the theoretical derivation step-by-step.
    Then, it uses a set of example parameters to calculate a numerical
    value for the critical speed, displaying the calculation process.
    """

    # --- Step 1: Theoretical Derivation ---
    print("Derivation of Critical Speed (v_crit) for an Oversteering Vehicle")
    print("="*70)
    print("This derivation uses the linear single-track (bicycle) model of a vehicle.")
    print("\nParameters:")
    print("  a   - distance from Center of Gravity (CG) to front axle (m)")
    print("  b   - distance from CG to rear axle (m)")
    print("  c_f - cornering stiffness of the front axle (N/rad)")
    print("  c_r - cornering stiffness of the rear axle (N/rad)")
    print("  m   - vehicle mass (kg)")
    print("  v   - vehicle forward speed (m/s)")
    print("-" * 70)

    print("\n1. Stability Condition:")
    print("The linear lateral dynamics of the vehicle become unstable when the following condition is met:")
    print("  det(A) <= 0")
    print("Where 'A' is the state matrix of the system. The boundary of stability (critical speed)")
    print("is found when det(A) = 0. This leads to the equation:")
    print("\n  (c_f + c_r)*(a^2*c_f + b^2*c_r) - (m*v^2 + a*c_f - b*c_r)*(a*c_f - b*c_r) = 0")

    print("\n2. Isolating v^2:")
    print("Rearranging the terms to solve for the speed (v), we get:")
    print("  m * v^2 * (a*c_f - b*c_r) = (c_f + c_r)*(a^2*c_f + b^2*c_r) - (a*c_f - b*c_r)^2")
    print("\nSimplifying the right-hand side yields a compact result:")
    print("  m * v^2 * (a*c_f - b*c_r) = c_f * c_r * (a+b)^2")
    
    print("\n3. Oversteering Requirement:")
    print("For a real solution for 'v' to exist, the vehicle must be oversteering.")
    print("This is defined by the condition: a*c_f > b*c_r")
    print("This makes the term (a*c_f - b*c_r) positive, ensuring v^2 is positive.")
    
    print("\n4. Final Formula for Critical Speed:")
    print("Solving for v, we get the critical speed, v_crit:")
    print("\n  v_crit = sqrt( (c_f * c_r * (a+b)^2) / (m * (a*c_f - b*c_r)) )")
    print("="*70)

    # --- Step 2: Numerical Calculation ---
    print("\nExample Calculation:")
    
    # Example parameters for an oversteering vehicle
    a = 1.1  # m
    b = 1.6  # m
    c_f = 100000  # N/rad
    c_r = 60000   # N/rad
    m = 1500      # kg
    
    # Check for oversteering condition
    oversteer_term = a * c_f - b * c_r
    if oversteer_term <= 0:
        print("\nThe provided parameters do not describe an oversteering vehicle.")
        print("A critical speed does not exist under these conditions.")
        return

    wheelbase = a + b
    numerator = c_f * c_r * (wheelbase**2)
    denominator = m * oversteer_term
    
    v_crit_sq = numerator / denominator
    v_crit = math.sqrt(v_crit_sq)

    print("\nUsing the following parameters:")
    print(f"  a   = {a} m")
    print(f"  b   = {b} m")
    print(f"  c_f = {c_f} N/rad")
    print(f"  c_r = {c_r} N/rad")
    print(f"  m   = {m} kg")
    
    print("\nThe calculation is:")
    # Per the instruction: "output each number in the final equation!"
    print(f"\nv_crit = sqrt( ({c_f} * {c_r} * ({a} + {b})^2) / ({m} * ({a}*{c_f} - {b}*{c_r})) )")
    print(f"v_crit = sqrt( ({c_f * c_r} * {wheelbase:.2f}^2) / ({m} * {oversteer_term}) )")
    print(f"v_crit = sqrt( {numerator} / {denominator} )")
    print(f"v_crit = sqrt( {v_crit_sq:.4f} )")
    print(f"\nFinal Result:")
    print(f"The critical speed is {v_crit:.4f} m/s.")
    print(f"This is equivalent to {v_crit * 3.6:.2f} km/h or {v_crit * 2.23694:.2f} mph.")

    # Return the final answer in the requested format
    global final_answer
    final_answer = v_crit

if __name__ == '__main__':
    final_answer = 0.0
    calculate_critical_speed()
    # The final answer format is handled outside of a typical script execution.
    # For this environment, we print it at the end.
    print(f'<<<{final_answer:.4f}>>>')
