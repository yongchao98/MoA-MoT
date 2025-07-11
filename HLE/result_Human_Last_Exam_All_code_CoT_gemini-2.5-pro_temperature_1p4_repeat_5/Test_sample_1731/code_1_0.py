import sympy
from sympy import init_printing, Integral, oo, pi, Symbol, Function, Eq, exp

def find_equilibrium_values():
    """
    This function uses symbolic mathematics to derive and display the equilibrium
    values for mean energy and entropy of a photon gas based on Bose-Einstein statistics.
    """
    # Initialize pretty printing for mathematical formulas
    init_printing(use_unicode=True)

    # --- Step 1: Define Symbols ---
    # Define physical constants and variables as symbols for the equations.
    k_B = Symbol('k_B')  # Boltzmann constant
    T = Symbol('T')      # Absolute temperature
    h = Symbol('h')      # Planck's constant
    c = Symbol('c')      # Speed of light in vacuum
    V = Symbol('V')      # Volume of the container
    E = Symbol('E')      # Total mean energy of the system
    S = Symbol('S')      # Total entropy of the system
    epsilon = Symbol('epsilon', positive=True) # Energy variable for integration

    print("--- Derivation of Equilibrium Values for a Photon Gas ---")
    print("\nBased on maximizing entropy for a Bose gas with non-conserved particle number,")
    print("we first find the equilibrium distribution of photons over energy states.")
    print("This yields the Bose-Einstein distribution for photons (where chemical potential is zero).")

    # --- Step 2: Calculate Equilibrium Mean Energy (E) ---
    print("\n1. Equilibrium Mean Energy (E)")
    print("The total energy E is found by integrating the energy of photons over the continuous")
    print("density of states g(epsilon) for a 3D volume V.")
    
    # Density of states for photons (including 2 polarization states)
    g_epsilon = (8 * pi * V) / (h**3 * c**3) * epsilon**2
    
    # Full integrand for energy
    integrand_E = epsilon * g_epsilon / (exp(epsilon / (k_B * T)) - 1)
    
    print("\nThe density of states g(epsilon) is given by:")
    print(sympy.pretty(Eq(Symbol('g(epsilon)'), g_epsilon)))

    print("\nThe total energy is the integral:")
    print(sympy.pretty(Eq(E, Integral(integrand_E, (epsilon, 0, oo)))))
    
    # The integral part is Integral(epsilon**3 / (exp(beta*epsilon) - 1)) with beta = 1/(k_B*T)
    # The standard form of this integral, Integral(x**3 / (e**x - 1) dx) from 0 to infinity,
    # is a known result from the study of the Riemann zeta function and is equal to pi**4 / 15.
    
    # Constructing the final expression for Energy
    energy_expr = (8 * pi**5 * V * (k_B*T)**4) / (15 * h**3 * c**3)
    equilibrium_energy_eq = Eq(E, energy_expr)
    
    print("\nEvaluating this integral gives the Stefan-Boltzmann law for the total energy of black-body radiation:")
    print("\n--- FINAL EQUATION FOR MEAN ENERGY (E) ---")
    print(sympy.pretty(equilibrium_energy_eq))
    print("-" * 45)

    # --- Step 3: Calculate Equilibrium Entropy (S) ---
    print("\n2. Equilibrium Entropy (S)")
    print("For a photon gas, there is a direct thermodynamic relationship between entropy S, energy E, pressure P, and volume V:")
    print("S = (E + PV) / T")
    print("The pressure of a photon gas is P = E / (3V). Substituting this gives:")
    print("S = (E + (E/3V)*V) / T = (4/3) * (E/T)")
    
    # Entropy in terms of Energy E and Temperature T
    entropy_expr_E = (4 / 3) * E / T
    equilibrium_entropy_eq = Eq(S, entropy_expr_E)

    print("\n--- FINAL EQUATION FOR ENTROPY (S) in terms of E ---")
    print(sympy.pretty(equilibrium_entropy_eq))
    print("-" * 55)

    # Substitute the expression for E to get S in terms of T
    entropy_expr_T = (4 / 3) * energy_expr / T
    equilibrium_entropy_eq_T = Eq(S, sympy.simplify(entropy_expr_T))

    print("\nBy substituting the expression for E, we can write S as a function of temperature T:")
    print("\n--- FINAL EQUATION FOR ENTROPY (S) in terms of T ---")
    print(sympy.pretty(equilibrium_entropy_eq_T))
    print("-" * 55)


if __name__ == '__main__':
    find_equilibrium_values()
