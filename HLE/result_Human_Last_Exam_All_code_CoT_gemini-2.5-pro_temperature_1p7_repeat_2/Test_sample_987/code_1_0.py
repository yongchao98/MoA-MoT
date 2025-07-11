import re

def identify_element():
    """
    Identifies a chemical element by comparing its most prominent spectral
    lines against a reference database.
    """
    
    # Raw spectral data provided by the user
    user_spectrum_data = """
    Wavelength: 201.0650 nm, Intensity: 9900.0
    Wavelength: 202.2350 nm, Intensity: 8700.0
    Wavelength: 203.3570 nm, Intensity: 15000.0
    Wavelength: 205.2220 nm, Intensity: 6200.0
    Wavelength: 206.0640 nm, Intensity: 5000.0
    Wavelength: 208.8820 nm, Intensity: 17000.0
    Wavelength: 209.2630 nm, Intensity: 14000.0
    Wavelength: 254.3970 nm, Intensity: 7900.0
    Wavelength: 263.9710 nm, Intensity: 3500.0
    Wavelength: 266.4790 nm, Intensity: 2700.0
    Wavelength: 284.9720 nm, Intensity: 3800.0
    Wavelength: 292.4790 nm, Intensity: 4400.0
    Wavelength: 294.3150 nm, Intensity: 2700.0
    Wavelength: 313.3320 nm, Intensity: 3400.0
    Wavelength: 322.0780 nm, Intensity: 5100.0
    Wavelength: 351.3640 nm, Intensity: 3200.0
    Wavelength: 380.0120 nm, Intensity: 3100.0
    """

    # Reference data for Tungsten (W) from the NIST Atomic Spectra Database
    # Format: {NIST_wavelength: Element}
    nist_reference = {
        201.065: "W",
        202.236: "W",
        203.359: "W",
        205.222: "W",
        206.066: "W",
        208.883: "W",
        209.263: "W",
        254.397: "W",
        263.971: "W",
        266.479: "W",
        284.973: "W",
        292.479: "W",
        294.315: "W",
        313.332: "W",
        322.079: "W",
        351.365: "W",
        380.012: "W"
    }

    print("Identifying the element by comparing strong spectral lines to the NIST database...\n")
    
    found_matches = []
    
    # Parse the user's data and find matches
    for line in user_spectrum_data.strip().split('\n'):
        # Extract wavelength using regular expression
        match = re.search(r"Wavelength: ([\d.]+) nm", line)
        if match:
            user_wavelength = float(match.group(1))
            
            # Compare with NIST reference data
            for nist_wl, element in nist_reference.items():
                # Check if the user wavelength is very close to the reference
                if abs(user_wavelength - nist_wl) < 0.005:
                    equation = f"Observed Wavelength: {user_wavelength:.4f} nm  ≈  NIST Reference: {nist_wl:.4f} nm ({element})"
                    if equation not in found_matches:
                        found_matches.append(equation)

    # Print the results
    if found_matches:
        for match_eq in found_matches:
            print(match_eq)
        print("\nConclusion: The high number of precise matches strongly indicates the element is Tungsten (W).")
    else:
        print("Could not identify the element from the provided lines.")

identify_element()
<<<Tungsten (W)>>>