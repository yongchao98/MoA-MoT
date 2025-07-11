import math

def derive_critical_speed():
    """
    Derives and calculates the critical speed for an oversteering vehicle
    based on the linear single-track model.
    """
    # Parameters of the linear single-track model
    # a: distance from CG to front axle (m)
    a = 1.1
    # b: distance from CG to rear axle (m)
    b = 1.5
    # c_f: cornering stiffness of the front axle (N/rad)
    c_f = 80000.0
    # c_r: cornering stiffness of the rear axle (N/rad)
    c_r = 50000.0
    # m: vehicle mass (kg)
    m = 1500.0
    # I: vehicle moment of inertia (kg*m^2) - not needed for this calculation but part of the model
    I = 3000.0

    print("Derivation of Critical Speed for an Oversteering Vehicle")
    print("="*60)
    print("The linear single-track model's stability is analyzed. The system becomes unstable")
    print("when the determinant of its state matrix goes to zero. Solving for speed 'v' at")
    print("this point gives the critical speed.")
    print("\nThe formula is: v_crit = sqrt( (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r)) )\n")
    print("Vehicle Parameters:")
    print(f"  a   = {a} m")
    print(f"  b   = {b} m")
    print(f"  c_f = {c_f} N/rad")
    print(f"  c_r = {c_r} N/rad")
    print(f"  m   = {m} kg\n")

    # Check for oversteering condition: a*c_f > b*c_r
    oversteer_check = a * c_f - b * c_r
    
    if oversteer_check <= 0:
        print("Stability Analysis:")
        print(f"The term (a*c_f - b*c_r) is {oversteer_check:.2f}.")
        print("The vehicle is not oversteering (it is understeering or neutral steering).")
        print("It is stable at all speeds and has no critical speed.")
        # We can't return the answer in this case, so we'll stop.
        return None
    else:
        print("Stability Analysis:")
        print(f"The term (a*c_f - b*c_r) = {a*c_f:.2f} - {b*c_r:.2f} = {oversteer_check:.2f} > 0.")
        print("The vehicle is oversteering, so a critical speed exists.\n")
        
        # Calculate the critical speed
        numerator = c_f * c_r * (a + b)**2
        denominator = m * oversteer_check
        v_crit_sq = numerator / denominator
        v_crit = math.sqrt(v_crit_sq)

        # Print the equation with numerical values substituted
        print("Final Equation with Numerical Values:")
        # The equation shows all the base parameters being substituted
        final_equation = (
            f"v_crit = sqrt( ({c_f} * {c_r} * ({a} + {b})^2) / ({m} * ({a} * {c_f} - {b} * {c_r})) )"
        )
        print(final_equation)
        
        # Print the final result
        print("\nResult:")
        print(f"The critical speed is: {v_crit:.4f} m/s")
        # Returning the value for the final answer block
        return v_crit

# Execute the function and capture the result
critical_speed_result = derive_critical_speed()
# The final answer must be in the required format.
if critical_speed_result is not None:
    # Use f-string to format the number without trailing zeros if it's an integer
    final_answer = f'{critical_speed_result:.4f}'
    # Remove unnecessary trailing zeros and the decimal point if it's an integer
    final_answer = final_answer.rstrip('0').rstrip('.')
    print(f'<<<{final_answer}>>>')
