import sys

def solve_mercury_tail_problem():
    """
    This script explains the visibility of Mercury's sodium tail
    when observed with a specific filter.
    """

    # Key parameters from the problem description
    filter_center_wavelength_nm = 486
    sodium_emission_wavelength_nm = 589

    # Explanation of the physics involved
    print("Step 1: Determine the emission source and color of Mercury's tail.")
    print(f"Mercury's tail is composed mainly of sodium atoms. These atoms emit light most strongly at the Sodium D-lines, which are at approximately {sodium_emission_wavelength_nm} nm. This wavelength corresponds to a yellow color.")
    print("-" * 30)

    print("Step 2: Determine what light the filter allows to pass.")
    print(f"The observation uses a narrow bandpass filter centered at {filter_center_wavelength_nm} nm. This filter is designed to let light around {filter_center_wavelength_nm} nm pass through, while blocking other wavelengths.")
    print("-" * 30)

    print("Step 3: Compare the filter's wavelength with the tail's emission wavelength.")
    print(f"The filter at {filter_center_wavelength_nm} nm is very far from the sodium emission at {sodium_emission_wavelength_nm} nm.")
    print("Therefore, the light that makes the sodium tail visible will be blocked by the filter.")
    print("-" * 30)

    print("Step 4: Conclude the final appearance.")
    print("Because the filter blocks the light from the sodium atoms, the sodium tail would not be visible in the image.")

solve_mercury_tail_problem()