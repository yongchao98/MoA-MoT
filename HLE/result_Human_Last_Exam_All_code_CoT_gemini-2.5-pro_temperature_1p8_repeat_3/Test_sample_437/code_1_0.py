import sympy
from sympy import symbols, exp, integrate, pi, oo, sqrt, sin, Mul, Pow

def solve_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a Slater-type orbital.
    """
    # Define symbols for spherical coordinates (r, theta, phi) and the orbital exponent (zeta)
    r, theta, phi = symbols('r theta phi')
    zeta = symbols('zeta', positive=True)

    # 1. Define the 1s Slater Orbital (unnormalized)
    # The functional form is exp(-zeta*r)
    phi_unnormalized = exp(-zeta * r)
    print(f"The unnormalized 1s Slater orbital is: phi(r) = {phi_unnormalized}")

    # 2. Normalize the orbital by solving <phi|phi> = 1
    # The normalization integral is Integral(N^2 * |phi|^2 * dtau) = 1
    # In spherical coordinates, dtau = r**2 * sin(theta) dr d(theta) d(phi)
    integrand_norm = phi_unnormalized**2 * r**2 * sin(theta)

    # The normalization integral is N^2 * Integral(integrand_norm)
    # We solve for N^2 = 1 / Integral(integrand_norm)
    norm_integral_val = integrate(integrand_norm, (r, 0, oo), (theta, 0, pi), (phi, 0, 2*pi))
    
    # N_squared is the square of the normalization constant
    N_squared = 1 / norm_integral_val
    print(f"\nThe square of the normalization constant, N^2, is: {N_squared}")
    
    # The normalization constant N
    N = sqrt(N_squared)
    print(f"The normalization constant, N, is: {N}")

    # The complete normalized 1s Slater orbital
    phi_normalized = N * phi_unnormalized
    print(f"\nThe full normalized 1s Slater orbital is: phi_1s(r) = {phi_normalized}")

    # 3. Evaluate the integral V = <phi_1s| 1/r |phi_1s>
    # The operator is 1/r
    operator = 1/r
    
    # The integrand is |phi_1s|^2 * (1/r) * dtau
    integrand_V11 = phi_normalized**2 * operator * r**2 * sin(theta)
    
    print("\nWe now evaluate the integral: <phi_1s| 1/r |phi_1s>")
    print("V = Integral( (phi_1s)^2 * (1/r) * r^2 * sin(theta) ) dr d(theta) d(phi)")
    
    # We can write V = N^2 * Integral(exp(-2*zeta*r) * r * sin(theta)) dr d(theta) d(phi)
    # Let's evaluate the components separately to build the final equation.

    # Angular integral: Integral(sin(theta) d(theta) d(phi))
    angular_integral = integrate(sin(theta), (theta, 0, pi), (phi, 0, 2*pi))

    # Radial integral: Integral(exp(-2*zeta*r) * r * dr)
    radial_integral = integrate(exp(-2*zeta*r) * r, (r, 0, oo))

    print("\n--- The Final Equation ---")
    print("The integral V is the product of three terms: (N^2) * (Radial Part) * (Angular Part)")
    
    # Constructing the equation with all the numbers/symbols
    term1 = N_squared
    term2 = radial_integral
    term3 = angular_integral
    
    print(f"Term 1 (N^2): {term1}")
    print(f"Term 2 (Radial Integral): {term2}")
    print(f"Term 3 (Angular Integral): {term3}")

    # The final equation showing all components before cancellation
    final_equation = Mul(term1, term2, term3, evaluate=False)
    print(f"\nFinal Equation: V = {term1} * ({term2}) * ({term3})")

    # 4. Calculate the final simplified result
    final_result = sympy.simplify(final_equation)

    print(f"\nAfter cancellation, the result simplifies to:")
    print(f"V = {final_result}")

if __name__ == '__main__':
    solve_integral()