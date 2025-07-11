import sympy as sp

def evaluate_slater_integral():
    """
    Symbolically evaluates the integral <phi_1s | 1/r | phi_1s>
    for a normalized 1s Slater-Type Orbital (STO).
    """
    
    # Define symbolic variables
    # r: radial distance
    # zeta: orbital exponent
    r = sp.symbols('r', real=True, positive=True)
    zeta = sp.symbols('zeta', real=True, positive=True, nonzero=True)
    
    # We work in spherical coordinates. The volume element is d(tau) = r^2 * sin(theta) dr d(theta) d(phi).
    # The integral of the angular part over the whole sphere (d(Omega)) is 4*pi.
    # integral(sin(theta) d(theta) from 0 to pi) * integral(d(phi) from 0 to 2*pi) = 2 * 2*pi = 4*pi
    angular_integral = 4 * sp.pi
    
    # --- 1. Find the Normalization Constant N ---
    # The normalization condition is Integral(|N*exp(-zeta*r)|^2) * d(tau) = 1.
    # N^2 * Integral(exp(-2*zeta*r) * r^2 dr) * Integral(d(Omega)) = 1.
    
    # The radial integral for normalization is:
    radial_integral_norm = sp.integrate(r**2 * sp.exp(-2 * zeta * r), (r, 0, sp.oo))
    
    # From N^2 * radial_integral_norm * angular_integral = 1, we can find N^2.
    N_squared = 1 / (radial_integral_norm * angular_integral)

    # --- 2. Evaluate the integral <phi_1s | 1/r | phi_1s> ---
    # This integral is Integral( (phi_1s)^2 * (1/r) * d(tau) )
    # = N^2 * Integral(exp(-2*zeta*r) * (1/r) * r^2 dr) * Integral(d(Omega))
    
    # The radial part of this expectation value integral is:
    radial_integral_exp = sp.integrate(r * sp.exp(-2 * zeta * r), (r, 0, sp.oo))

    # The final result is the product of the three components:
    final_result = N_squared * radial_integral_exp * angular_integral

    # --- 3. Print the results step-by-step ---
    print("Evaluation of the integral <phi_1s | 1/r | phi_1s> for a 1s Slater orbital.")
    print("-" * 70)
    print("The integral is calculated as the product of three terms:")
    print("1. The squared normalization constant (N^2)")
    print("2. The radial part of the expectation value integral")
    print("3. The angular part of the integral (solid angle)")
    print("-" * 70)
    
    print(f"Squared Normalization Constant (N^2): {N_squared}")
    print(f"Radial Integral (Integral r * exp(-2*zeta*r) dr): {radial_integral_exp}")
    print(f"Angular Integral (Integral d(Omega)): {angular_integral}")

    print("\n" + "-" * 70)
    print("Final Equation:")
    print(f"Result = (N^2) * (Radial Integral) * (Angular Integral)")
    # The line below shows each component in the final multiplication
    print(f"({N_squared}) * ({radial_integral_exp}) * ({angular_integral}) = {final_result}")
    print("-" * 70)

if __name__ == '__main__':
    evaluate_slater_integral()
