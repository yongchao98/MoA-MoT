import math

def calculate_critical_speed():
    """
    Calculates and prints the critical speed for an oversteering vehicle
    based on the linear single-track model.
    """
    # Parameters of the linear single-track model
    # Distance from CG to front axle (m)
    a = 1.5
    # Distance from CG to rear axle (m)
    b = 1.0
    # Cornering stiffness of the front axle (N/rad)
    c_f = 90000.0
    # Cornering stiffness of the rear axle (N/rad)
    c_r = 60000.0
    # Vehicle mass (kg)
    m = 1500.0
    # Note: Moment of inertia 'I' is part of the full dynamic model but
    # cancels out or is not present in the final stability condition (det(A)=0) equation.

    # Check for oversteering condition (a*c_f > b*c_r)
    # This is required for a real critical speed to exist.
    if a * c_f <= b * c_r:
        print("The provided parameters do not describe an oversteering vehicle.")
        print(f"Condition for oversteer (a*c_f > b*c_r) is not met:")
        print(f"{a} * {c_f} = {a * c_f}")
        print(f"{b} * {c_r} = {b * c_r}")
        return

    # Calculate the components of the formula
    wheelbase_sq = (a + b)**2
    numerator = c_f * c_r * wheelbase_sq
    denominator = m * (a * c_f - b * c_r)

    # Calculate critical speed squared and then the critical speed
    v_crit_sq = numerator / denominator
    v_crit = math.sqrt(v_crit_sq)

    # Output the final equation with numerical values
    print("Derivation of the critical speed (v_crit) for an oversteering vehicle.")
    print("\nFormula: v_crit = sqrt( (c_f * c_r * (a+b)^2) / (m * (a*c_f - b*c_r)) )")
    
    print("\nSubstituting the given parameters:")
    print(f"a   = {a} m")
    print(f"b   = {b} m")
    print(f"c_f = {c_f} N/rad")
    print(f"c_r = {c_r} N/rad")
    print(f"m   = {m} kg")
    
    # Building the equation string with numbers
    # To maintain clarity, let's calculate the intermediate terms for the printout
    term_a_cf = a * c_f
    term_b_cr = b * c_r
    term_a_plus_b_sq = (a + b)**2
    
    equation = (
        f"v_crit = sqrt( ({c_f} * {c_r} * ({a} + {b})^2) / ({m} * ({a} * {c_f} - {b} * {c_r})) )\n"
        f"       = sqrt( ({c_f * c_r} * {term_a_plus_b_sq:.2f}) / ({m} * ({term_a_cf} - {term_b_cr})) )\n"
        f"       = sqrt( {numerator:.2e} / ({m} * {term_a_cf - term_b_cr}) )\n"
        f"       = sqrt( {numerator:.2e} / {denominator:.2e} )\n"
        f"       = sqrt( {v_crit_sq:.2f} )\n"
    )

    print("\nFinal Equation:")
    print(equation)

    print(f"The critical speed is {v_crit:.2f} m/s.")
    # You can uncomment the line below to convert it to km/h
    # print(f"Which is equal to {v_crit * 3.6:.2f} km/h.")
    
    # Final answer as per requested format
    global final_answer
    final_answer = v_crit

# Execute the function
calculate_critical_speed()
# The required format is <<<number>>>, so we print the final_answer variable.
print(f'<<<{final_answer:.2f}>>>')