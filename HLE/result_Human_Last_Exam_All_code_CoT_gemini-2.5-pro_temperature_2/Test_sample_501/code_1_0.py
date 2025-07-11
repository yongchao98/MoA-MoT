import math

def calculate_isolated_polymer_force(n, l, x, E_0):
    """
    Calculates the restorative force of a thermally isolated freely jointed polymer chain.

    The force law is derived using the microcanonical ensemble, where entropy is conserved.
    f(x) = -[3*x*E_0 / (n * l^2 * (n-1))] * exp[3*x^2 / (2 * n * l^2 * (n-1))]

    Args:
        n (int): The number of monomers in the chain (assumed to be large).
        l (float): The length of a single strut.
        x (float): The end-to-end separation of the polymer.
        E_0 (float): The kinetic energy of the polymer at zero extension (x=0).

    Returns:
        float: The restorative force exerted by the polymer.
    """
    # The derivation is valid for x << n*l
    if x >= n * l:
        print("Warning: The extension 'x' is not small compared to the total contour")
        print("length 'n*l'. The Gaussian approximation may not be accurate.")

    # Calculate the terms in the force equation
    prefactor_numerator = 3 * x * E_0
    prefactor_denominator = n * (l**2) * (n - 1)
    
    exponent_numerator = 3 * x**2
    exponent_denominator = 2 * n * (l**2) * (n - 1)

    prefactor = prefactor_numerator / prefactor_denominator
    exponent = exponent_numerator / exponent_denominator

    force = -prefactor * math.exp(exponent)
    
    # --- Outputting the results as requested ---
    print("The force law for the thermally isolated polymer is:")
    print("f(x) = - (3 * x * E(0)) / (n * l^2 * (n-1)) * exp( (3 * x^2) / (2 * n * l^2 * (n-1)) )\n")

    print("For the given values:")
    print(f"  n (number of monomers)        = {n}")
    print(f"  l (strut length)              = {l}")
    print(f"  x (end-to-end separation)     = {x}")
    print(f"  E(0) (initial kinetic energy) = {E_0:.4g}\n")

    print("The components of the equation are calculated as follows:")
    
    print("\n1. Prefactor: -(3 * x * E(0)) / (n * l^2 * (n-1))")
    print(f"   Numerator (3 * x * E(0)):      3 * {x} * {E_0:.4g} = {prefactor_numerator:.4g}")
    print(f"   Denominator (n * l^2 * (n-1)): {n} * {l}^2 * {n-1} = {prefactor_denominator:.4g}")
    print(f"   Resulting Prefactor:            -({prefactor_numerator:.4g} / {prefactor_denominator:.4g}) = {-prefactor:.4g}")

    print("\n2. Exponential Term: exp( (3 * x^2) / (2 * n * l^2 * (n-1)) )")
    print(f"   Numerator (3 * x^2):         3 * {x}^2 = {exponent_numerator:.4g}")
    print(f"   Denominator (2*n*l^2*(n-1)): 2 * {n} * {l}^2 * {n-1} = {exponent_denominator:.4g}")
    print(f"   Exponent value:              {exponent_numerator:.4g} / {exponent_denominator:.4g} = {exponent:.4g}")
    print(f"   exp(exponent):                 exp({exponent:.4g}) = {math.exp(exponent):.4g}")
    
    print("\nFinal Force Calculation:")
    print(f"f(x) = (Prefactor) * exp(exponent)")
    print(f"f(x) = ({-prefactor:.4g}) * ({math.exp(exponent):.4g})")
    print(f"f({x}) = {force:.4g}")

    return force

# Example usage with physically plausible values
# n: number of monomers
# l: length in meters (e.g., Kuhn length of dsDNA is ~100 nm)
# x: extension in meters
# E_0: Initial energy, let's use the average energy at 300K for ν degrees of freedom
# E_0 = (nu/2) * k_B * T ≈ (n-1) * k_B * T
k_B = 1.380649e-23  # Boltzmann constant in J/K
T = 300            # Temperature in Kelvin
n_val = 1000
l_val = 100e-9     # 100 nm
x_val = 1e-6       # 1 µm extension (1% of total length)
E_0_val = (n_val-1) * k_B * T

_ = calculate_isolated_polymer_force(n=n_val, l=l_val, x=x_val, E_0=E_0_val)