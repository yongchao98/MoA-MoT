import sympy as sp

def display_absorption_equation(case):
    """
    Constructs and prints the equation for the absorption cross-section
    for a given physical scenario.

    Args:
        case (str): 'a' for isolated molecules, 'b' for interacting molecules.
    """
    # Define the symbolic variables for the equations
    sigma = sp.Symbol('sigma(omega_0)') # Absorption cross-section
    omega_0 = sp.Symbol('omega_0')      # Laser center frequency
    tau = sp.Symbol('tau')              # Gaussian pulse duration
    hbar = sp.Symbol('hbar')            # Reduced Planck constant
    c = sp.Symbol('c')                  # Speed of light in vacuum
    n_refr = sp.Symbol('n')             # Refractive index of the medium
    epsilon_0 = sp.Symbol('epsilon_0')  # Vacuum permittivity
    pi = sp.pi

    print("-" * 70)
    
    if case.lower() == 'a':
        # --- Case a: Isolated Molecules ---
        omega_fi = sp.Symbol('omega_fi')   # Single-molecule transition frequency
        mu_fi_sq = sp.Symbol('|mu_fi|^2')  # Squared transition dipole moment of a molecule

        # Assemble the equation for case (a)
        equation = (2 * sp.sqrt(pi) * tau * omega_fi * mu_fi_sq) / (hbar * c * n_refr * epsilon_0) * sp.exp(-(omega_fi - omega_0)**2 * tau**2)
        
        # Print explanation and the final formatted equation
        print("Case a) Non-Interacting Molecules\n")
        print("Description: Molecules are treated as independent absorbers. The transition is between single-molecule electronic states.")
        print("\nComponents of the equation:")
        print(f"  sigma(omega_0): The absorption cross-section as a function of the laser frequency.")
        print(f"  omega_0:       The central frequency of the Gaussian laser pulse.")
        print(f"  omega_fi:      The transition frequency of a single molecule, (E_final - E_initial)/hbar.")
        print(f"  |mu_fi|^2:      The squared transition dipole moment of the molecule.")
        print(f"  tau:           The duration of the Gaussian laser pulse.")
        print(f"  n:             The refractive index of the medium.")
        print(f"  c, hbar, epsilon_0: Fundamental constants.")
        
        print("\nThe equation for the absorption cross-section is:")
        sp.pprint(sp.Eq(sigma, equation), use_unicode=True)

    elif case.lower() == 'b':
        # --- Case b: Interacting Molecules (Frenkel Exciton Model) ---
        omega_exc = sp.Symbol('omega_exc')  # Exciton transition frequency
        mu_ge_sq = sp.Symbol('|mu_ge|^2')   # Squared transition dipole moment of a single molecule
        N = sp.Symbol('N')                  # Number of molecules in the chain
        J = sp.Symbol('J')                  # Near-neighbor coupling energy
        E_mol = sp.Symbol('E_mol')          # Single-molecule excitation energy

        # Exciton dipole moment is coherently enhanced
        mu_exc_sq = N * mu_ge_sq
        
        # Assemble the equation for case (b)
        equation = (2 * sp.sqrt(pi) * tau * omega_exc * mu_exc_sq) / (hbar * c * n_refr * epsilon_0) * sp.exp(-(omega_exc - omega_0)**2 * tau**2)
        
        # Print explanation and the final formatted equation
        print("Case b) Interacting Molecules (Frenkel Exciton Model)\n")
        print("Description: Near-neighbor interactions create delocalized excitons. Absorption is dominated by the k=0 bright exciton state.")
        print("\nComponents of the equation:")
        print(f"  sigma(omega_0): The absorption cross-section of the exciton peak.")
        print(f"  omega_0:       The central frequency of the Gaussian laser pulse.")
        print(f"  omega_exc:     The transition frequency of the k=0 bright exciton, given by ({E_mol} + 2*{J})/{hbar}.")
        print(f"  N:             The number of interacting molecules in the chain.")
        print(f"  |mu_ge|^2:      The squared transition dipole moment for a single molecule.")
        print(f"  Note: The total dipole strength is N*|mu_ge|^2 (superradiance).")
        print(f"  tau:           The duration of the Gaussian laser pulse.")
        print(f"  n:             The refractive index of the medium.")
        print(f"  c, hbar, epsilon_0: Fundamental constants.")

        print("\nThe equation for the absorption cross-section is:")
        sp.pprint(sp.Eq(sigma, equation), use_unicode=True)
    else:
        print(f"Error: Case '{case}' is not recognized. Please use 'a' or 'b'.")
        
    print("-" * 70)

if __name__ == '__main__':
    # --- Execute and Print Results ---
    # Call the function for Case a
    display_absorption_equation('a')
    
    # Call the function for Case b
    display_absorption_equation('b')