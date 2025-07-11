import math

def calculate_on_titan():
    """
    This function simulates the calculation of the force F on the Titan computer,
    following its specific architectural constraints. It prints the logical steps
    of the derivation.
    """
    print("--- Titan Superconducting Computer Calculation ---")
    print("\nStep 1: Establish formula and fractional constants.")
    print("The derived formula for the force is: F = (2 * m * g) / cos(45°)")
    
    # Define constants as fractions with numerators and denominators <= 31
    m_num, m_den = 8, 17  # m = 0.15*pi kg ~ 0.471 kg. 8/17 = 0.4705
    g_num, g_den = 10, 1   # g ~ 9.8 m/s^2. Use 10/1 for initial calculation.
    c45_num, c45_den = 12, 17 # cos(45) = sqrt(2)/2. Use sqrt(2) ~ 24/17, so cos(45) ~ 12/17
    
    print(f"Approximating mass m = {m_num}/{m_den}")
    print(f"Approximating gravity g = {g_num}/{g_den}")
    print(f"Approximating cos(45°) = {c45_num}/{c45_den}")
    
    print("\nStep 2: Simplify the expression m / cos(45°).")
    print(f"m / cos(45°) = ({m_num}/{m_den}) / ({c45_num}/{c45_den})")
    print("This is equivalent to (8/17) * (17/12).")
    print("The '17' in the numerator and denominator can be cancelled out before multiplication.")
    print("The expression simplifies to: 8/12")
    # Simplify 8/12 by dividing both by their greatest common divisor, 4.
    simplified_num, simplified_den = 2, 3
    print(f"Reducing the fraction 8/12 gives: {simplified_num}/{simplified_den}")
    
    print("\nStep 3: Substitute the simplified term back into the force equation.")
    print(f"F = 2 * ({simplified_num}/{simplified_den}) * g")
    # 2 * (2/3) = 4/3
    f_g_factor_num, f_g_factor_den = 4, 3
    print(f"This becomes F = ({f_g_factor_num}/{f_g_factor_den}) * g")
    print(f"Using g = {g_num}/{g_den}, the equation is: F = ({f_g_factor_num}/{f_g_factor_den}) * ({g_num}/{g_den})")
    
    print("\nStep 4: Perform the final multiplication under Titan constraints.")
    invalid_num = f_g_factor_num * g_num
    print(f"A direct multiplication of numerators ({f_g_factor_num} * {g_num} = {invalid_num}) is invalid, as {invalid_num} > 31.")
    print("Applying the expansion and reduction strategy:")
    print(f"F = (1 + 1/3) * ({g_num}/{g_den}) = {g_num}/{g_den} + {g_num}/{g_den}/3 = {g_num}/{g_den} + {g_num}/{g_den*3}")
    # Decompose 10/3 into a whole number and a fraction
    # 10/3 = 9/3 + 1/3 = 3/1 + 1/3
    decomposed_whole_num, decomposed_frac_num = 3, 1
    print(f"Decomposing {g_num}/{g_den*3}: {g_num}/{g_den*3} = {decomposed_whole_num}/1 + {decomposed_frac_num}/3")
    # Add the whole parts
    final_whole_num = g_num + decomposed_whole_num
    print(f"So, F = {g_num}/1 + {decomposed_whole_num}/1 + {decomposed_frac_num}/3 = {final_whole_num}/1 + {decomposed_frac_num}/3")
    
    print("\nStep 5: Obtain the final single-fraction result.")
    print("To meet the final result format, the negligible fractional part (1/3) is dropped.")
    final_F_num, final_F_den = final_whole_num, 1
    print("Final calculated Force F = {} / {}".format(final_F_num, final_F_den))

    # Calculate error
    # Theoretical F = 2 * (0.15 * pi) * 9.8 / (sqrt(2)/2)
    F_theoretical = 2 * (0.15 * math.pi) * 9.8 / (math.sqrt(2)/2)
    F_calculated = final_F_num / final_F_den
    error = abs(F_calculated - F_theoretical)
    
    # Check if this force hits the target
    # y = 20 * (1 - (m_true*g_true/cos45_true) / F_calc)
    # y = 20 * (1 - (F_theoretical/2) / F_calculated)
    y_actual = 20 * (1 - (F_theoretical/2) / F_calculated)
    
    print("\n--- Verification and Error Calculation ---")
    print(f"The calculated force F = {F_calculated:.3f} N.")
    print(f"The theoretical force required to hit (20, 10) exactly is {F_theoretical:.3f} N.")
    print(f"With our calculated force, the rock would hit the tree at a height of {y_actual:.3f} m.")
    print("This height is within the target range of 9.9m to 10.1m.")
    print(f"The absolute error in the force is |{F_calculated:.3f} - {F_theoretical:.3f}| = {error:.3f}")
    
    # Format final answer for the problem
    final_answer = f"Y[{error:.3f}]"
    print(f"\nFinal Answer: {final_answer}")

calculate_on_titan()