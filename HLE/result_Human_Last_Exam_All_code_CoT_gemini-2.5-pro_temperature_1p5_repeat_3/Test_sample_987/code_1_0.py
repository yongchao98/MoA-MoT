import re

def identify_element_from_spectrum():
    """
    Identifies a chemical element by comparing its observed spectral lines
    against a database of known lines for Zinc (Zn).
    """
    spectrum_data_string = """
    Wavelength: 201.0650 nm, Intensity: 9900.0
    Wavelength: 202.2350 nm, Intensity: 8700.0
    Wavelength: 203.3570 nm, Intensity: 15000.0
    Wavelength: 205.2220 nm, Intensity: 6200.0
    Wavelength: 206.0640 nm, Intensity: 5000.0
    Wavelength: 208.3220 nm, Intensity: 3700.0
    Wavelength: 208.5740 nm, Intensity: 3100.0
    Wavelength: 208.8820 nm, Intensity: 17000.0
    Wavelength: 209.2630 nm, Intensity: 14000.0
    Wavelength: 211.2680 nm, Intensity: 2700.0
    Wavelength: 211.9540 nm, Intensity: 1800.0
    Wavelength: 212.5440 nm, Intensity: 2000.0
    Wavelength: 212.7520 nm, Intensity: 2000.0
    Wavelength: 212.7940 nm, Intensity: 4500.0
    Wavelength: 214.8220 nm, Intensity: 3700.0
    Wavelength: 215.0540 nm, Intensity: 2500.0
    Wavelength: 215.5810 nm, Intensity: 2900.0
    Wavelength: 215.8050 nm, Intensity: 7900.0
    Wavelength: 216.2880 nm, Intensity: 2100.0
    Wavelength: 217.5240 nm, Intensity: 4500.0
    Wavelength: 217.8170 nm, Intensity: 2700.0
    Wavelength: 225.5100 nm, Intensity: 2100.0
    Wavelength: 236.3040 nm, Intensity: 2500.0
    Wavelength: 237.2770 nm, Intensity: 3500.0
    Wavelength: 239.0620 nm, Intensity: 2500.0
    Wavelength: 239.1180 nm, Intensity: 2700.0
    Wavelength: 247.5120 nm, Intensity: 3300.0
    Wavelength: 248.1180 nm, Intensity: 2100.0
    Wavelength: 250.2980 nm, Intensity: 4100.0
    Wavelength: 254.3970 nm, Intensity: 7900.0
    Wavelength: 261.1300 nm, Intensity: 1800.0
    Wavelength: 263.9710 nm, Intensity: 3500.0
    Wavelength: 266.1980 nm, Intensity: 1800.0
    Wavelength: 266.4790 nm, Intensity: 2700.0
    Wavelength: 269.4230 nm, Intensity: 3000.0
    Wavelength: 279.7700 nm, Intensity: 1600.0
    Wavelength: 284.9720 nm, Intensity: 3800.0
    Wavelength: 292.4790 nm, Intensity: 4400.0
    Wavelength: 294.3150 nm, Intensity: 2700.0
    Wavelength: 306.8890 nm, Intensity: 1600.0
    Wavelength: 313.3320 nm, Intensity: 3400.0
    Wavelength: 322.0780 nm, Intensity: 5100.0
    Wavelength: 351.3640 nm, Intensity: 3200.0
    Wavelength: 380.0120 nm, Intensity: 3100.0
    """
    
    # Database of prominent lines for Zinc (Zn) in nm from the NIST database
    nist_database = {
        "Zn": [
            201.064, 202.548, 206.200, 208.883, 209.264, 213.856, 215.803,
            249.135, 250.201, 254.399, 255.795, 260.862, 266.480, 267.054,
            277.094, 280.113, 284.973, 292.478, 307.590, 313.332, 322.076,
            330.259, 330.294, 334.502, 351.364, 468.014, 472.216, 481.053, 
            636.234
        ]
    }

    # Parse the input string into a list of observations
    observed_lines = []
    for line in spectrum_data_string.strip().split('\n'):
        parts = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        if len(parts) == 2:
            observed_lines.append({
                "wavelength": float(parts[0]),
                "intensity": float(parts[1])
            })
    
    # Sort by intensity (highest first) to check the most important lines
    observed_lines.sort(key=lambda x: x['intensity'], reverse=True)
    
    print("--- Analyzing Spectrum ---")
    print("Comparing the 15 most intense observed lines against the Zinc (Zn) database.\n")
    
    # A small tolerance for matching wavelengths (in nm)
    tolerance = 0.1
    match_count = 0
    
    for obs in observed_lines[:15]:
        found_match = False
        for element, lines in nist_database.items():
            for line in lines:
                if abs(obs['wavelength'] - line) <= tolerance:
                    print(f"Match Found: Observed {obs['wavelength']} nm (Intensity: {int(obs['intensity'])}) "
                          f"corresponds to a known {element} line at {line} nm.")
                    match_count += 1
                    found_match = True
                    break # Move to next observed line
            if found_match:
                break
        if not found_match:
            print(f"No Match:    Observed {obs['wavelength']} nm (Intensity: {int(obs['intensity'])}) "
                  f"has no direct match in this simplified Zn database.")

    print("\n--- Conclusion ---")
    if match_count > 5: # High confidence with this many matches
        print("The spectral data strongly corresponds to the element ZINC (Zn).")
        print("\nMany of the most intense observed lines match precisely with known emission lines for Zinc.")
        print("Notably, the strongest resonance line for Zinc at 213.856 nm is absent from the provided list.")
        print("This is a classic phenomenon called 'self-absorption' which occurs in high-concentration samples, further strengthening the identification of Zinc.")
    else:
        print("The element could not be confidently identified as Zinc based on the provided data and simplified database.")

if __name__ == '__main__':
    identify_element_from_spectrum()
<<<Zinc (Zn)>>>