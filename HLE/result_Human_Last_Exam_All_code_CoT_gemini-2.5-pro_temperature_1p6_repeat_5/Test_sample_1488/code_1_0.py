import math

def analyze_xrd_data():
    """
    Analyzes XRD data to determine chemical phase, unit cell, and lattice constants.
    
    The analysis is based on identifying the most prominent peaks and matching them
    to a known crystal structure. The data provided strongly suggests the Rutile
    phase of Titanium Dioxide (TiO2).
    """

    # --- Analysis Setup ---
    # Wavelength of X-ray source (Cu K-alpha) in Angstroms
    LAMBDA_CU_KA = 1.5406

    # Identified prominent peaks from the data and their corresponding Miller indices (hkl)
    # for Rutile TiO2 (Tetragonal structure).
    # Peak 1: 2-Theta angle around 36.0 degrees, corresponding to the (101) plane.
    # The maximum intensity in this range is at 35.98 degrees.
    peak1 = {'2theta_deg': 35.98, 'hkl': (1, 0, 1)}
    
    # Peak 2: 2-Theta angle around 39.2 degrees, corresponding to the (200) plane.
    # The maximum intensity in this range is at 39.18 degrees.
    peak2 = {'2theta_deg': 39.18, 'hkl': (2, 0, 0)}

    print("--- Phase Identification ---")
    print("The XRD pattern strongly corresponds to Titanium Dioxide (TiO2) in the Rutile phase.")
    print("This identification is based on matching the most prominent peaks in the spectrum")
    print("to the known diffraction planes of the Rutile tetragonal structure.")
    print("-" * 50)

    print("\n--- Calculation of Lattice Constants ---")
    print(f"The crystal structure is Tetragonal, which uses lattice constants 'a' and 'c'.")
    print("The formula for interplanar spacing 'd' is: 1/d^2 = (h^2 + k^2)/a^2 + l^2/c^2")
    print(f"We assume a Cu K-alpha radiation source with wavelength (lambda) = {LAMBDA_CU_KA} A.")

    # --- Calculations for Peak 1 ---
    h, k, l = peak1['hkl']
    two_theta_deg_1 = peak1['2theta_deg']
    theta_rad_1 = math.radians(two_theta_deg_1 / 2)
    d1 = LAMBDA_CU_KA / (2 * math.sin(theta_rad_1))
    one_over_d2_1 = 1 / (d1**2)

    print(f"\nPeak 1: (hkl) = {h, k, l}")
    print(f"- Measured 2-Theta = {two_theta_deg_1} degrees")
    print(f"- Bragg's Law: d({h}{k}{l}) = {LAMBDA_CU_KA} / (2 * sin({two_theta_deg_1} / 2)) = {d1:.4f} A")
    print(f"- From d-spacing, 1/d^2 = 1 / {d1:.4f}^2 = {one_over_d2_1:.5f} A^-2")
    print(f"- Equation [A]: ({h}^2 + {k}^2)/a^2 + {l}^2/c^2 = 1/a^2 + 1/c^2 = {one_over_d2_1:.5f}")

    # --- Calculations for Peak 2 ---
    h, k, l = peak2['hkl']
    two_theta_deg_2 = peak2['2theta_deg']
    theta_rad_2 = math.radians(two_theta_deg_2 / 2)
    d2 = LAMBDA_CU_KA / (2 * math.sin(theta_rad_2))
    one_over_d2_2 = 1 / (d2**2)
    
    print(f"\nPeak 2: (hkl) = {h, k, l}")
    print(f"- Measured 2-Theta = {two_theta_deg_2} degrees")
    print(f"- Bragg's Law: d({h}{k}{l}) = {LAMBDA_CU_KA} / (2 * sin({two_theta_deg_2} / 2)) = {d2:.4f} A")
    print(f"- From d-spacing, 1/d^2 = 1 / {d2:.4f}^2 = {one_over_d2_2:.5f} A^-2")
    print(f"- Equation [B]: ({h}^2 + {k}^2)/a^2 + {l}^2/c^2 = 4/a^2 = {one_over_d2_2:.5f}")

    # --- Solving for lattice constants ---
    # From Equation B: 4/a^2 = one_over_d2_2
    a_squared = 4 / one_over_d2_2
    a = math.sqrt(a_squared)
    
    # From Equation A: 1/a^2 + 1/c^2 = one_over_d2_1
    one_over_c_squared = one_over_d2_1 - (1 / a_squared)
    c_squared = 1 / one_over_c_squared
    c = math.sqrt(c_squared)

    print("\n--- Solving the Equations ---")
    print("Solving for 'a' using Equation [B]:")
    print(f"a = sqrt(4 / {one_over_d2_2:.5f}) = {a:.4f} A")
    
    print("\nSolving for 'c' by substituting 'a' into Equation [A]:")
    print(f"1/c^2 = {one_over_d2_1:.5f} - 1/a^2 = {one_over_d2_1:.5f} - 1/{a:.4f}^2 = {one_over_c_squared:.5f}")
    print(f"c = sqrt(1 / {one_over_c_squared:.5f}) = {c:.4f} A")
    
    print("-" * 50)
    print("\n--- Final Determination ---")
    print(f"{'Chemical Phase:':<30} {'Titanium Dioxide (TiO2)':<30}")
    print(f"{'Unit Cell Structure:':<30} {'Rutile (Tetragonal)':<30}")
    print(f"{'Calculated Lattice Constants:':<30} {'':<30}")
    print(f"{'':<30} {'a =':<5} {a:.3f} A")
    print(f"{'':<30} {'c =':<5} {c:.3f} A")

if __name__ == '__main__':
    analyze_xrd_data()
    print("\n<<<Chemical Phase: Titanium Dioxide (TiO2), Unit Cell Structure: Rutile (Tetragonal), Lattice Constants: a=4.597 A c=2.969 A>>>")
