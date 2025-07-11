import xraylib
import math

def calculate_transmission():
    """
    Calculates the transmission of characteristic X-rays through a Beryllium window
    for a given list of elements.
    """
    # --- Parameters of the EDX setup ---
    # Beryllium window properties
    window_material_symbol = 'Be'
    window_material_name = 'Beryllium'
    # Density of Beryllium in g/cm^3
    rho_be = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(window_material_symbol))
    # Thickness of the window in cm (100 µm = 0.01 cm)
    t_be = 100.0 * 1e-4

    # --- Elements to check from the answer choices ---
    # We evaluate them in order of atomic number (lightest first)
    elements_to_check = ['Na', 'Mg', 'Si', 'Ca']

    print(f"Calculating X-ray transmission through a {t_be*1e4} µm {window_material_name} window.\n")

    for symbol in elements_to_check:
        try:
            # Get element's atomic number
            Z = xraylib.SymbolToAtomicNumber(symbol)

            # Get the energy of the K-alpha line in keV
            # Note: X-ray energies are in keV in xraylib, but CS functions expect keV.
            E_keV = xraylib.LineEnergy(Z, xraylib.KA_LINE)

            # Get the mass attenuation coefficient (mu) of Be at this energy (in cm^2/g)
            mu_be = xraylib.CS_Total_CP(window_material_symbol, E_keV)

            # Calculate the exponent term in the Beer-Lambert law
            exponent = mu_be * rho_be * t_be

            # Calculate the transmission T = exp(-exponent)
            transmission = math.exp(-exponent)

            # Print the results for each step of the calculation
            print(f"Element: {symbol} (Z={Z})")
            print(f"  Kα Energy (E): {E_keV:.4f} keV")
            print(f"  Be Window: Thickness (t) = {t_be:.4f} cm, Density (ρ) = {rho_be:.4f} g/cm³")
            print(f"  Be Attenuation at E (μ): {mu_be:.4f} cm²/g")
            print(f"  Final Equation: Transmission = exp(-(μ * ρ * t))")
            print(f"  Equation Values: T = exp(-({mu_be:.4f} * {rho_be:.4f} * {t_be:.4f}))")
            print(f"  Exponent (μ * ρ * t): {exponent:.4f}")
            print(f"  Transmission: {transmission * 100:.6f} %")
            print("-" * 40)

        except ValueError as e:
            print(f"Could not perform calculation for {symbol}: {e}")
            print("-" * 40)

if __name__ == '__main__':
    calculate_transmission()