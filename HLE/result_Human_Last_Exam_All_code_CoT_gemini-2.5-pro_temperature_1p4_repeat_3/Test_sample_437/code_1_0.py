import sympy as sp

def solve_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital
    and prints the step-by-step derivation.
    """
    # Define symbols
    # zeta is the orbital exponent, a positive real number.
    # r is the radial distance, a positive real number.
    zeta = sp.Symbol('zeta', positive=True)
    r = sp.Symbol('r', positive=True)
    theta, phi = sp.symbols('theta phi')

    print("Evaluating the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital.")
    print("-" * 70)

    # 1. The normalized 1s Slater-Type Orbital (STO)
    # The normalization constant N is found by solving Integral[|phi_1s|^2 dV] = 1.
    # For phi_1s = N * exp(-zeta*r), this gives N^2 = zeta**3/pi.
    N_sq = zeta**3 / sp.pi
    
    print("The integral is I = Integral[ phi_1s^* * (1/r) * phi_1s * dV ] over all space.")
    print("For a real 1s orbital phi_1s(r) = N * exp(-zeta*r), this becomes:")
    print("I = Integral[ N^2 * exp(-2*zeta*r) * (1/r) * dV ]")
    print(f"The squared normalization constant N^2 = {N_sq}")
    
    # 2. Separate into angular and radial parts in spherical coordinates
    # dV = r^2 * sin(theta) * dr * dtheta * dphi
    print("\nIn spherical coordinates, the integral separates into two parts:")
    
    # 3. Evaluate the angular part
    # The integral of sin(theta) over theta (0 to pi) and phi (0 to 2*pi).
    angular_integral = sp.integrate(sp.sin(theta), (phi, 0, 2 * sp.pi), (theta, 0, sp.pi))
    print("\n1. Angular Part: Integral[ sin(theta) dtheta dphi ]")
    print(f"   Value of the angular integral = {angular_integral}")

    # 4. Evaluate the radial part
    # The radial part is Integral[ N^2 * exp(-2*zeta*r) * (1/r) * r^2 dr ] from 0 to infinity.
    # This simplifies to N^2 * Integral[ r * exp(-2*zeta*r) dr ].
    print("\n2. Radial Part: Integral[ N^2 * r * exp(-2*zeta*r) dr ]")
    
    # Evaluate the core radial integral first
    core_radial_integral = sp.integrate(r * sp.exp(-2 * zeta * r), (r, 0, sp.oo))
    print(f"   The integral of r*exp(-2*zeta*r) from r=0 to infinity is: {core_radial_integral}")

    # Combine with N^2
    full_radial_integral = N_sq * core_radial_integral
    print(f"   The full radial part of the integral is (N^2) * ({core_radial_integral}) = {full_radial_integral}")

    # 5. Combine the parts to get the final result
    total_integral = full_radial_integral * angular_integral
    print("\n3. Total Integral = (Radial Part) * (Angular Part)")
    print(f"   Equation: 1 * <phi_1s| 1/r |phi_1s> = ({full_radial_integral}) * ({angular_integral})")
    print(f"   Symbolic Result = {total_integral}")

    print("-" * 70)
    print(f"The general evaluation of the integral <phi_1s| 1/r |phi_1s> is zeta.")

    # Substitute a numerical value for zeta for a specific case.
    # For a hydrogen atom, the orbital exponent zeta is 1.
    print("\nFor the specific case of a hydrogen atom, we set the orbital exponent zeta = 1.")
    zeta_val = 1
    numerical_result = total_integral.subs(zeta, zeta_val)
    print(f"For zeta = {zeta_val}, the value of the integral is {numerical_result}.")

if __name__ == '__main__':
    solve_integral()