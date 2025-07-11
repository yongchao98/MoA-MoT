import sympy as sp

def evaluate_1s_integral():
    """
    This function symbolically evaluates the integral <phi|1/r|phi> for a 1s Slater orbital
    and prints the step-by-step derivation.
    """
    # Define the necessary mathematical symbols.
    # zeta is the orbital exponent, and r is the radial distance.
    # Both are positive real numbers.
    zeta = sp.Symbol('zeta', positive=True)
    r = sp.Symbol('r', positive=True)
    N = sp.Symbol('N', positive=True)  # The normalization constant

    print("Step-by-step evaluation of <phi_1s| 1/r |phi_1s>\n")
    print("--------------------------------------------------")
    print("The 1s Slater-type orbital (STO) is given by: phi_1s = N * exp(-zeta*r)")
    print("where N is the normalization constant and zeta is the orbital exponent.\n")

    # Part 1: Determine the normalization constant N
    # The normalization condition is that the integral of |phi_1s|^2 over all space is 1.
    # For a spherically symmetric function, the volume element dV in spherical coordinates
    # integrates over the angular parts to 4*pi, so dV = 4*pi*r^2 dr.
    # Equation: Integral_from_0_to_inf [ (N * exp(-zeta*r))^2 * 4*pi*r^2 ] dr = 1
    print("1. Normalization of the 1s orbital")
    print("   The condition is Integral(|phi_1s|^2 dV) = 1")
    print("   => Integral_0^inf [ (N * exp(-zeta*r))^2 * 4*pi*r^2 ] dr = 1")
    print("   => 4 * pi * N^2 * Integral_0^inf [ r^2 * exp(-2*zeta*r) ] dr = 1\n")

    # Evaluate the radial integral required for normalization.
    norm_integral = sp.integrate(r**2 * sp.exp(-2 * zeta * r), (r, 0, sp.oo))
    print(f"   Let's evaluate the integral part: Integral_0^inf [ r^2 * exp(-2*zeta*r) ] dr")
    print(f"   Using standard integrals, the result is: {norm_integral}\n")

    # Substitute the integral result back into the normalization equation to solve for N^2.
    N_squared_expr = 1 / (4 * sp.pi * norm_integral)
    print("   Substituting this back into the normalization equation:")
    print(f"   4 * pi * N^2 * ({norm_integral}) = 1")
    # Simplify the expression for N^2
    N_squared = sp.simplify(N_squared_expr)
    print(f"   Solving for N^2 gives: N^2 = 1 / (4 * pi * {norm_integral}) = {N_squared}\n")

    # Part 2: Evaluate the expectation value integral <phi_1s| 1/r |phi_1s>
    # Integral I = Integral(phi_1s^* * (1/r) * phi_1s dV)
    # Since phi_1s is a real function, phi_1s^* = phi_1s.
    # I = Integral_0^inf [ (N * exp(-zeta*r))^2 * (1/r) * 4*pi*r^2 ] dr
    # I = 4 * pi * N^2 * Integral_0^inf [ r * exp(-2*zeta*r) ] dr
    print("2. Evaluation of the integral <phi_1s| 1/r |phi_1s>")
    print("   The integral is I = Integral( phi_1s * (1/r) * phi_1s * dV )")
    print("   => I = Integral_0^inf [ (N*exp(-zeta*r))^2 * (1/r) * 4*pi*r^2 ] dr")
    print("   which simplifies to:")
    print("   => I = 4 * pi * N^2 * Integral_0^inf [ r * exp(-2*zeta*r) ] dr\n")

    # Substitute the derived value of N^2 into the expression for I.
    I_coeff = 4 * sp.pi * N_squared
    print(f"   We found N^2 = {N_squared}. Substituting this into the equation for I:")
    # Using the requirement to output each number (or symbol)
    print(f"   I = ({4}) * ({sp.pi}) * ({N_squared}) * Integral_0^inf [ r * exp(-2*zeta*r) ] dr")
    print(f"   I = {sp.simplify(I_coeff)} * Integral_0^inf [ r * exp(-2*zeta*r) ] dr\n")

    # Evaluate the final radial integral.
    potential_integral = sp.integrate(r * sp.exp(-2 * zeta * r), (r, 0, sp.oo))
    print(f"   Now, we evaluate the remaining integral: Integral_0^inf [ r * exp(-2*zeta*r) ] dr")
    print(f"   The result of this integral is: {potential_integral}\n")

    # Combine all parts to get the final result.
    final_result = sp.simplify(I_coeff * potential_integral)

    # Part 3: Display the Final Result
    print("3. Final Result")
    print("   Combining all the calculated parts:")
    # Display the final equation with its components before the final simplification.
    print(f"   <phi_1s| 1/r |phi_1s> = ({sp.simplify(I_coeff)}) * ({potential_integral})")
    print(f"   which computes to:")
    print(f"   I = {final_result}\n")
    print("--------------------------------------------------")
    print(f"The final evaluated value of the integral is: {final_result}")

if __name__ == '__main__':
    evaluate_1s_integral()