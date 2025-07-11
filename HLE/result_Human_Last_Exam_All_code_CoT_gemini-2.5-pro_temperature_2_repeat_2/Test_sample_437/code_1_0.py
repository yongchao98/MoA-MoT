def evaluate_slater_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a Slater-type orbital
    and shows the calculation step-by-step.
    """
    print("This script evaluates the expectation value of the 1/r operator for a 1s Slater orbital.")
    print("The integral is: <phi_1s | 1/r | phi_1s> = integral(phi_1s * (1/r) * phi_1s * d_tau)\n")
    print("---------------------------------------------------------------------")

    print("Step 1: Define the components")
    print("The normalized 1s Slater-type orbital is: phi_1s = (zeta^3 / pi)^(1/2) * exp(-zeta * r)")
    print("The operator is: 1/r")
    print("The volume element in spherical coordinates is: d_tau = r^2 * sin(theta) * dr * d(theta) * d(phi)\n")

    print("The full integral expression is:")
    print("I = integral( (zeta^3/pi) * exp(-2*zeta*r) * (1/r) * r^2 * sin(theta) ) dr d(theta) d(phi)")
    print("\nThis expression can be separated into parts:")
    print("I = (zeta^3 / pi) * [integral from 0 to inf of r*exp(-2*zeta*r) dr] * [integral from 0 to pi of sin(theta) d(theta)] * [integral from 0 to 2*pi of d(phi)]")
    print("---------------------------------------------------------------------")

    print("Step 2: Evaluate the separated integrals")

    # Angular Part
    angular_phi_result_str = "2*pi"
    angular_theta_result_str = "2"
    total_angular_result_str = "4*pi"
    print(f"The angular part of the integral evaluates to:")
    print(f"  - Integral d(phi) from 0 to 2*pi = {angular_phi_result_str}")
    print(f"  - Integral sin(theta) d(theta) from 0 to pi = {angular_theta_result_str}")
    print(f"  - Total angular contribution = (2*pi) * 2 = {total_angular_result_str}\n")

    # Radial Part
    radial_result_str = "1 / (4 * zeta^2)"
    print(f"The radial part is the integral of r * exp(-2*zeta*r) dr from 0 to infinity.")
    print(f"Using the standard integral formula integral(x^n * exp(-ax) dx) = n! / a^(n+1) with n=1 and a=2*zeta:")
    print(f"  - The result is 1! / (2*zeta)^2 = {radial_result_str}\n")

    # Normalization Constant squared
    norm_const_squared_str = "zeta^3 / pi"
    print(f"The squared normalization constant from the two wavefunctions is: {norm_const_squared_str}")
    print("---------------------------------------------------------------------")

    print("Step 3: Combine all parts to get the final result")
    print("Result = (Squared Normalization) * (Angular Part) * (Radial Part)")

    final_equation = f"Result = ({norm_const_squared_str}) * ({total_angular_result_str}) * ({radial_result_str})"
    print("Final Equation: " + final_equation)
    print("Result = zeta\n")
    
    print("The numerical constants and powers in the final combined equation are:")
    # The equation is: (zeta^3 / pi) * (4 * pi) * (1 / (4 * zeta^2))
    print(f"Power of zeta in normalization factor: 3")
    print(f"Numerical factor from angular integrals: 4")
    print(f"Numerator from radial integral: 1")
    print(f"Power of zeta from radial integral: 2")

if __name__ == '__main__':
    evaluate_slater_integral()
