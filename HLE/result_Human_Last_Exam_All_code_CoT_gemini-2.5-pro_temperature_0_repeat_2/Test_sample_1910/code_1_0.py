import math
from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator

def find_second_peak_q_space():
    """
    Calculates the Q-space position of the second major diffraction peak
    for NaMgH3 perovskite.
    """
    # Step 1: Define the crystal structure of NaMgH3 (Pnma, No. 62)
    # Lattice parameters and atomic positions are taken from experimental data
    # (e.g., Materials Project entry mp-29913, based on literature).
    try:
        structure = Structure.from_spacegroup(
            "Pnma",
            lattice=[[5.530, 0, 0], [0, 7.808, 0], [0, 0, 5.514]],
            species=["Na", "Mg", "H", "H"],
            coords=[
                [0.016, 0.25, 0.981],  # Na on Wyckoff site 4c
                [0, 0, 0.5],           # Mg on Wyckoff site 4b
                [0.578, 0.25, 0.081],  # H1 on Wyckoff site 4c
                [0.298, 0.041, 0.298]   # H2 on Wyckoff site 8d
            ]
        )
    except ImportError as e:
        print(f"Error: Pymatgen library not found. Please install it using 'pip install pymatgen'")
        print(f"Details: {e}")
        return
    except Exception as e:
        print(f"An error occurred during structure creation: {e}")
        return

    # Step 2: Simulate the powder diffraction pattern.
    # The wavelength is given as 0.2952 Å.
    wavelength = 0.2952
    calculator = XRDCalculator(wavelength=wavelength)
    pattern = calculator.get_pattern(structure)

    # Step 3: Combine pattern data and sort by intensity to find major peaks.
    # The pattern contains 2-theta, intensity, Miller indices (hkl), and d-spacing.
    peaks = []
    for i in range(len(pattern.y)):
        intensity = pattern.y[i]
        # Group degenerate peaks, but treat them as one for sorting
        if intensity > 0.1: # Filter out very weak peaks
            peaks.append({
                'intensity': intensity,
                'd_hkl': pattern.d_hkls[i],
                'hkl': pattern.hkls[i][0]['hkl'] # Get Miller indices of the first degenerate peak
            })

    # Sort peaks by intensity in descending order
    sorted_peaks = sorted(peaks, key=lambda p: p['intensity'], reverse=True)

    # The second major peak is the second element in the sorted list.
    second_peak = sorted_peaks[1]
    d_hkl = second_peak['d_hkl']
    hkl = second_peak['hkl']

    # Step 4: Calculate the Q-space position for the second peak.
    q_value = 2 * math.pi / d_hkl

    # Print the results and the calculation steps
    print(f"The second major diffraction peak corresponds to the Miller indices (hkl) = {hkl}.")
    print(f"The d-spacing for this peak is {d_hkl:.4f} Å.")
    print("\nThe Q-space position is calculated using the formula: Q = 2 * pi / d")
    print(f"Q = 2 * {math.pi:.6f} / {d_hkl:.4f}")
    print(f"Q = {q_value:.4f} 1/Å")
    
    # Return the final answer in the specified format
    print(f"\n<<<{q_value:.4f}>>>")

if __name__ == '__main__':
    find_second_peak_q_space()