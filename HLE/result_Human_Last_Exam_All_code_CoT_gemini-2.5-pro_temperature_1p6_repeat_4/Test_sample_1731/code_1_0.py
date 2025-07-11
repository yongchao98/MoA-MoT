import numpy as np

def calculate_photon_gas_properties():
    """
    Calculates and prints the equilibrium properties of a photon gas.
    
    This script details the formulas for the equilibrium distribution, mean energy density,
    and entropy density of a photon gas, derived from the principles of statistical mechanics.
    It then calculates the numerical constants in the energy and entropy equations.
    """

    # --- Physical Constants ---
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck constant in J*s
    c = 299792458     # Speed of light in m/s

    # --- Introduction and Equilibrium Distribution ---
    print("For a gas of light quanta (photons), the non-conservation of particle number sets the chemical potential to zero.")
    print("The equilibrium distribution for the average occupation number n(epsilon) of a state with energy epsilon is the Planck Distribution:")
    print("Equation [1]: n(epsilon) = 1 / (exp(epsilon / (k_B * T)) - 1)\n")

    # --- Equilibrium Mean Energy Calculation ---
    print("The equilibrium mean energy density (u) is given by the Stefan-Boltzmann Law:")
    print("Equation [2]: u = a * T^4")
    
    # Calculate the radiation constant 'a'
    # a = (8 * pi^5 * k_B^4) / (15 * h^3 * c^3)
    a = (8 * np.pi**5 * k_B**4) / (15 * h**3 * c**3)
    
    print(f"The radiation constant 'a' is calculated to be {a:.4e} J m^-3 K^-4.")
    print("Thus, the final equation for mean energy density is:")
    print(f"u = {a:.4e} * T^4  [J/m^3]\n")

    # --- Equilibrium Entropy Calculation ---
    print("The equilibrium entropy density (s) is derived from the energy density:")
    print("Equation [3]: s = (4/3) * a * T^3")

    # Calculate the coefficient for the entropy density equation
    entropy_coeff = (4/3) * a

    print(f"The coefficient for entropy density (4/3)*a is calculated to be {entropy_coeff:.4e} J m^-3 K^-4.")
    print("Thus, the final equation for entropy density is:")
    print(f"s = {entropy_coeff:.4e} * T^3  [J/m^3 K]")

if __name__ == '__main__':
    calculate_photon_gas_properties()