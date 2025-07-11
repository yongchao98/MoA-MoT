import numpy as np
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator

def find_second_peak_q_space():
    """
    Finds the Q-space position for the second major diffraction peak of NaMgH3.
    """
    # Step 1: Get the crystal structure from the Materials Project database.
    # NaMgH3 at room temperature has material ID "mp-22923".
    material_id = "mp-22923"
    print(f"Fetching crystal structure for NaMgH3 ({material_id}) from the Materials Project...")
    try:
        with MPRester() as mpr:
            structure = mpr.get_structure_by_material_id(material_id)
        print("Successfully fetched structure.")
        print(f"Formula: {structure.formula}")
        print(f"Lattice parameters (a, b, c): {structure.lattice.abc}")
    except Exception as e:
        print(f"Could not fetch structure. Error: {e}")
        print("Using fallback data for NaMgH3 (orthorhombic, Pnma).")
        # Fallback in case the API fails.
        # This data is from Materials Project for mp-22923.
        from pymatgen.core import Structure
        lattice = [[5.6171, 0, 0], [0, 7.9419, 0], [0, 0, 5.5684]]
        species = ["Na", "Na", "Na", "Na", "Mg", "Mg", "Mg", "Mg", 
                   "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"]
        coords = [[0.007, 0.25, 0.463], [0.507, 0.75, 0.963], [0.993, 0.75, 0.537], [0.493, 0.25, 0.037],
                  [0.0, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.0],
                  [0.575, 0.25, 0.165], [0.075, 0.75, 0.665], [0.425, 0.75, 0.835], [0.925, 0.25, 0.335],
                  [0.280, 0.038, 0.720], [0.780, 0.538, 0.220], [0.220, 0.962, 0.280], [0.720, 0.462, 0.780],
                  [0.720, 0.962, 0.220], [0.220, 0.538, 0.780], [0.780, 0.038, 0.280], [0.280, 0.462, 0.720]]
        structure = Structure(lattice, species, coords)


    # Step 2: Calculate the diffraction pattern.
    # The wavelength value affects the 2-theta angles but not the d-spacings or relative intensities.
    wavelength = 0.2952  # Angstrom
    calculator = XRDCalculator(wavelength=wavelength)
    pattern = calculator.get_pattern(structure)

    # Step 3: Identify the second major peak by sorting all peaks by intensity.
    # We combine the peak data into a list of dictionaries to sort it easily.
    all_peaks = []
    for i in range(len(pattern.x)):
        # For peaks with the same 2-theta (overlapping), pymatgen lists multiple hkls.
        # We take the first hkl listed as representative.
        hkl = pattern.hkls[i][0]['hkl']
        d_hkl = pattern.d_hkls[i]
        intensity = pattern.y[i]
        all_peaks.append({'hkl': hkl, 'd_hkl': d_hkl, 'intensity': intensity})

    # Sort by intensity in descending order
    sorted_peaks = sorted(all_peaks, key=lambda p: p['intensity'], reverse=True)

    first_peak = sorted_peaks[0]
    second_peak = sorted_peaks[1]

    print("\n--- Top Diffraction Peaks (by intensity) ---")
    print(f"1st Major Peak: hkl={first_peak['hkl']}, d-spacing={first_peak['d_hkl']:.4f} Å, Relative Intensity={first_peak['intensity']:.2f}")
    print(f"2nd Major Peak: hkl={second_peak['hkl']}, d-spacing={second_peak['d_hkl']:.4f} Å, Relative Intensity={second_peak['intensity']:.2f}")
    
    # Step 4: Calculate Q-space position for the second peak.
    d_second_peak = second_peak['d_hkl']
    q_second_peak = 2 * np.pi / d_second_peak
    
    print("\n--- Final Calculation ---")
    print(f"The second major peak has an interplanar spacing d = {d_second_peak:.4f} Å.")
    print("The Q-space position is calculated as Q = 2 * pi / d.")
    print(f"Q = 2 * {np.pi:.5f} / {d_second_peak:.4f} Å = {q_second_peak:.4f} Å⁻¹")
    print("\n-------------------------------------------")
    print(f"The second major diffraction peak is located at Q = {q_second_peak:.4f} 1/Angstrom.")

find_second_peak_q_space()
>>> 2.7437