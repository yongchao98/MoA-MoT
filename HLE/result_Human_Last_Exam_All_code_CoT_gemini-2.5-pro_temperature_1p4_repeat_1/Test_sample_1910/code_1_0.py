# The script requires the 'pymatgen' library.
# You can install it by running: pip install pymatgen

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.analysis.diffraction.xrd import XRDCalculator

def find_diffraction_peak_q_space():
    """
    This function calculates the Q-space position of the second major diffraction
    peak for NaMgH3 perovskite.
    """
    # Step 1: Define the crystal structure of NaMgH3 at room temperature.
    # Data from Bodd et al., Inorg. Chem. 2008, 47, 4391-4395.
    a = 5.6748  # Angstrom
    b = 7.9575  # Angstrom
    c = 5.5516  # Angstrom
    lattice = Lattice.orthorhombic(a, b, c)

    # Define the species and coordinates of the asymmetric unit.
    species = ["Na", "Mg", "H", "H"]
    # Coordinates for Na(4c), Mg(4b), H1(4c), H2(8d) sites.
    coords = [
        [0.0351, 0.25, 0.9702],    # Na
        [0.0, 0.0, 0.5],            # Mg
        [0.4691, 0.25, 0.1082],    # H1
        [0.7029, 0.0401, 0.7225]    # H2
    ]

    # Create the pymatgen Structure object using the space group "Pnma".
    structure = Structure.from_spacegroup("Pnma", lattice, species, coords)

    # Step 2: Simulate the X-ray diffraction pattern.
    # The wavelength is provided but not directly needed for Q if we have d.
    # Pymatgen's XRDCalculator uses it to determine angle-dependent factors
    # for a more accurate intensity calculation.
    wavelength_A = 0.2952
    xrd_calculator = XRDCalculator(wavelength=wavelength_A)
    pattern = xrd_calculator.get_pattern(structure, scaled=True, two_theta_range=None)

    # Step 3: Identify the second major peak by intensity.
    # We combine the calculated intensities, d-spacings, and hkl indices.
    # The pattern.y attribute holds intensities, and pattern.d_hkls holds d-spacings.
    peaks_data = sorted(zip(pattern.y, pattern.d_hkls, pattern.hkls), key=lambda x: x[0], reverse=True)

    # The most intense peak is the first element, the second most intense is the second.
    if len(peaks_data) < 2:
        print("Error: Less than two diffraction peaks were found.")
        return

    second_major_peak = peaks_data[1]

    # Step 4: Extract d-spacing and calculate Q.
    d_spacing_2nd = second_major_peak[1]
    hkl_2nd = second_major_peak[2][0]['hkl'] # Miller indices for the peak

    # The Q-space position is calculated using Q = 2 * pi / d.
    Q_2nd = 2 * np.pi / d_spacing_2nd

    print("The second major diffraction peak for NaMgH3 corresponds to the Miller indices (hkl):", hkl_2nd)
    print("The d-spacing for this peak is calculated to be {:.5f} Å.".format(d_spacing_2nd))
    print("\nThe Q-space position (Q) is calculated using the formula Q = 2 * π / d:")
    # Final output showing the numbers in the equation
    print(f"Q = 2 * {np.pi:.5f} / {d_spacing_2nd:.5f}")
    print(f"Q = {Q_2nd:.5f} Å⁻¹")

if __name__ == "__main__":
    find_diffraction_peak_q_space()