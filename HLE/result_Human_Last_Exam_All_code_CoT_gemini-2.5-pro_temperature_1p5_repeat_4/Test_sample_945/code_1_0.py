import math

def calculate_critical_speed():
    """
    Derives and calculates the critical speed for an oversteering vehicle
    based on the linear single-track model.
    """
    # --- Parameters of the linear single-track model ---
    # a: distance from CG to front axle (m)
    a = 1.2
    # b: distance from CG to rear axle (m)
    b = 1.6
    # c_f: cornering stiffness of the front axle (N/rad)
    c_f = 120000.0
    # c_r: cornering stiffness of the rear axle (N/rad)
    c_r = 80000.0
    # m: vehicle mass (kg)
    m = 1500.0
    # I: vehicle moment of inertia about the vertical axis (kg*m^2)
    # This parameter affects damping but not the critical speed itself.
    I = 3000.0

    print("Vehicle Parameters:")
    print(f"  a (CG to front axle distance): {a} m")
    print(f"  b (CG to rear axle distance): {b} m")
    print(f"  c_f (Front axle cornering stiffness): {c_f} N/rad")
    print(f"  c_r (Rear axle cornering stiffness): {c_r} N/rad")
    print(f"  m (Vehicle mass): {m} kg")
    print(f"  I (Vehicle moment of inertia): {I} kg*m^2\n")

    # --- Check for Oversteer Condition ---
    # The critical speed exists only for oversteering vehicles.
    # Oversteer condition: a*c_f > b*c_r
    oversteer_term_front = a * c_f
    oversteer_term_rear = b * c_r

    if oversteer_term_front <= oversteer_term_rear:
        print("Stability Analysis:")
        print(f"The condition for oversteer (a*c_f > b*c_r) is not met.")
        print(f"  a*c_f = {oversteer_term_front:.2f}")
        print(f"  b*c_r = {oversteer_term_rear:.2f}")
        print("The vehicle is understeering or neutral steer, which is stable at all speeds.")
        print("Therefore, a critical speed for divergence does not exist.")
        return

    # --- Derivation and Calculation ---
    print("Derivation of Critical Speed (v_crit):")
    print("The critical speed is the velocity at which the vehicle's lateral dynamics become unstable.")
    print("For an oversteering vehicle, this occurs when v > v_crit.")
    print("The formula is derived by finding the speed at which the constant term of the system's characteristic equation becomes zero.\n")
    
    # Numerator of the term inside the square root
    numerator = c_f * c_r * (a + b)**2
    # Denominator of the term inside the square root
    denominator = m * (oversteer_term_front - oversteer_term_rear)

    v_crit = math.sqrt(numerator / denominator)

    # --- Print the final equation with numerical values ---
    print("Final Equation with substituted values:")
    equation = f"v_crit = sqrt( (c_f * c_r * (a + b)^2) / (m * (a*c_f - b*c_r)) )"
    print(equation)
    
    equation_with_values = (
        f"v_crit = sqrt( ({c_r} * {c_f} * ({a} + {b})^2) / "
        f"({m} * ({a}*{c_f} - {b}*{c_r})) )"
    )
    print(equation_with_values)

    equation_calculated_values = (
        f"v_crit = sqrt( ({numerator:.2e}) / "
        f"({m} * ({oversteer_term_front:.2f} - {oversteer_term_rear:.2f})) )"
    )
    print(equation_calculated_values)
    
    equation_final_values = (
        f"v_crit = sqrt( {numerator:.4e} / {denominator:.4e} )"
    )
    print(equation_final_values)


    print("\n--- Result ---")
    print(f"The critical speed is: {v_crit:.2f} m/s")
    print(f"This is equivalent to {v_crit * 3.6:.2f} km/h.")

    # Final answer for the system
    global final_answer
    final_answer = v_crit

if __name__ == '__main__':
    final_answer = None
    calculate_critical_speed()
    if final_answer is not None:
        print(f"<<<{final_answer}>>>")
