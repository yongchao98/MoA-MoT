import sympy
from sympy import symbols, integrate, exp, pi, oo, sin, sqrt, pretty_print

def evaluate_integral():
    """
    Symbolically evaluates the integral <φ_1s| 1/r |φ_1s> for a 1s Slater-type orbital.
    The code will print the components of the calculation as per the user request.
    """
    # Define symbolic variables for the integration
    r, theta, phi = symbols('r theta phi', real=True, positive=True)
    zeta = symbols('zeta', real=True, positive=True)

    print("Evaluating the integral: <φ_1s| 1/r |φ_1s>\n")

    # 1. Define the normalized 1s Slater Orbital (φ_1s)
    # The normalization constant N is sqrt(zeta**3 / pi)
    N = sqrt(zeta**3 / pi)
    
    print(f"The normalized 1s Slater orbital is φ_1s = N * exp(-ζ*r), where N = {N}")
    
    # The full integral is ∫ (φ_1s)² * (1/r) * dV
    # dV = r² * sin(θ) dr dθ dφ
    # The integrand simplifies to N² * r * exp(-2ζr) * sin(θ)
    # We can separate this into three parts for the calculation.

    # 2. Calculate each part of the integral
    
    # Part 1: The square of the normalization constant
    N_squared = N**2
    
    # Part 2: The angular integral: ∫dφ from 0 to 2π and ∫sin(θ)dθ from 0 to π
    angular_integral = integrate(sin(theta), (phi, 0, 2 * pi), (theta, 0, pi))

    # Part 3: The radial integral: ∫ r * exp(-2ζr) dr from 0 to ∞
    radial_integrand = r * exp(-2 * zeta * r)
    radial_integral = integrate(radial_integrand, (r, 0, oo))

    # 3. Combine the parts to get the final result
    # Result = N² * (Angular Integral) * (Radial Integral)
    final_result = N_squared * angular_integral * radial_integral
    
    # 4. Print the final equation with each component's value
    print("\nThe final result is the product of three parts: (N²) * (Angular Integral) * (Radial Integral)")
    print("The final equation is:")
    print(f"Result = ({N_squared}) * ({angular_integral}) * ({radial_integral})")

    # 5. Print the simplified final answer
    print(f"\nSimplifying the expression gives:")
    print(f"Result = {final_result}")


evaluate_integral()