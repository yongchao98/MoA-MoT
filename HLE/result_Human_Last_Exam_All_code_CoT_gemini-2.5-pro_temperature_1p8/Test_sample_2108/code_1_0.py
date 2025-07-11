import sympy
from sympy import sin, cos, pi, integrate, Symbol, sqrt, latex

def solve_radiation_ratio():
    """
    Calculates the maximum achievable ratio of bidirectional conical power
    to line intensity for the given rotating sphere problem.
    """
    # Define mathematical symbols
    I_0 = Symbol('I_0', positive=True) # Base intensity constant
    theta = Symbol('theta')            # Polar angle
    phi = Symbol('phi')                # Azimuthal angle

    # The time-averaged intensity pattern from the rotating magnetic dipole is I(theta) = I_0 * (1 + cos(theta)^2)
    intensity_function = I_0 * (1 + cos(theta)**2)

    # To maximize the ratio P_cone / I_line, we must maximize the numerator and minimize the denominator.

    # 1. Minimize the denominator (I_line).
    # The minimum intensity occurs when theta = pi/2.
    min_line_intensity = intensity_function.subs(theta, pi/2)

    # 2. Maximize the numerator (P_cone).
    # This is done by aligning the bidirectional cone (half-angle pi/4) with the rotation axis (theta=0).
    # The power is integrated over two opposite cones. By symmetry, this is twice the integral over one cone.
    # dOmega = sin(theta) d(theta) d(phi)
    power_in_one_cone = integrate(intensity_function * sin(theta), (phi, 0, 2 * pi), (theta, 0, pi / 4))
    max_cone_power = 2 * power_in_one_cone

    # 3. Calculate the maximum achievable ratio.
    # The I_0 factor cancels out.
    max_ratio = max_cone_power / min_line_intensity

    # --- Output the results step-by-step ---
    print("Step 1: Define the intensity function and the ratio.")
    print(f"Intensity I(theta) = I_0 * (1 + cos(theta)^2)")
    print(f"Ratio to maximize: R = P_cone / I_line")
    print("-" * 30)

    print("Step 2: Calculate the minimum line intensity (denominator).")
    print("This occurs at theta = pi/2, perpendicular to the rotation axis.")
    print(f"I_line_min = I_0 * (1 + cos(pi/2)^2)")
    # We use '.doit()' to evaluate the symbolic expression.
    print(f"Denominator Value = {min_line_intensity.doit()}")
    print("-" * 30)
    
    print("Step 3: Calculate the maximum cone power (numerator).")
    print("This occurs by aligning the pi/4 half-angle cone with the rotation axis.")
    # Show the components of the calculation.
    integral_val_theta = integrate((1 + cos(theta)**2) * sin(theta), (theta, 0, pi / 4))
    integral_val_phi = integrate(1, (phi, 0, 2*pi))
    
    print(f"P_cone_max = 2 * I_0 * [Integral d(phi) from 0 to 2*pi] * [Integral (1+cos^2(t))sin(t) dt from 0 to pi/4]")
    print(f"Integral over phi = {integral_val_phi}")
    print(f"Integral over theta = {integral_val_theta}")
    print(f"Numerator Value = 2 * I_0 * ({integral_val_phi}) * ({integral_val_theta})")
    print(f"Numerator Value = {max_cone_power.doit()}")
    print("-" * 30)

    print("Step 4: Final maximum ratio calculation.")
    print("Ratio = (P_cone_max) / (I_line_min)")
    
    # Using symbolic representation for clarity in the equation.
    numerator_expr_str = latex(max_cone_power).replace("I_{0}", "")
    denominator_expr_str = latex(min_line_intensity).replace("I_{0}", "")

    # Output the final equation with all numbers and its evaluated result.
    print("\n--- Final Equation ---")
    print(f"Ratio = [ {latex(max_cone_power)} ] / [ {latex(min_line_intensity)} ]")
    print(f"Ratio = I_0 * ({numerator_expr_str}) / (I_0 * {denominator_expr_str})")
    print(f"Ratio = {max_ratio}")
    print(f"Numerically, the ratio is approximately: {max_ratio.evalf()}")


solve_radiation_ratio()
<<<6.38838634335436>>>