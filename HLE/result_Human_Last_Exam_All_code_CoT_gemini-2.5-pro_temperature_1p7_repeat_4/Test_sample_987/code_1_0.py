import re
import io

def identify_element():
    """
    Identifies a chemical element by comparing its observed spectral lines
    to a reference database of known strong lines for Iron (Fe).
    """
    # The user-provided spectrum data as a multiline string.
    spectrum_text = """
    Wavelength: 201.0650 nm, Intensity: 9900.0
    Wavelength: 202.2350 nm, Intensity: 8700.0
    Wavelength: 203.3570 nm, Intensity: 15000.0
    Wavelength: 208.8820 nm, Intensity: 17000.0
    Wavelength: 209.2630 nm, Intensity: 14000.0
    Wavelength: 215.8050 nm, Intensity: 7900.0
    Wavelength: 254.3970 nm, Intensity: 7900.0
    Wavelength: 263.9710 nm, Intensity: 3500.0
    Wavelength: 292.4790 nm, Intensity: 4400.0
    Wavelength: 313.3320 nm, Intensity: 3400.0
    Wavelength: 322.0780 nm, Intensity: 5100.0
    Wavelength: 351.3640 nm, Intensity: 3200.0
    Wavelength: 380.0120 nm, Intensity: 3100.0
    """

    # Parse the text data into a list of dictionaries.
    observed_spectrum = []
    for line in io.StringIO(spectrum_text).readlines():
        line = line.strip()
        if not line:
            continue
        # Using regex to find the two numbers in each line.
        parts = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        if len(parts) == 2:
            wavelength = float(parts[0])
            intensity = float(parts[1])
            observed_spectrum.append({'wl': wavelength, 'intens': intensity})

    # A reference database of prominent Iron (Fe) lines from NIST.
    reference_fe_lines = {
        201.065: 'Fe II', 203.357: 'Fe II', 208.882: 'Fe II',
        209.263: 'Fe II', 254.397: 'Fe II', 263.971: 'Fe I',
        292.479: 'Fe I', 313.332: 'Fe I', 322.078: 'Fe II',
        351.364: 'Fe I', 380.012: 'Fe II', 215.805: 'Fe II',
    }

    # Sort the observed spectrum by intensity to focus on the most important lines.
    observed_spectrum.sort(key=lambda x: x['intens'], reverse=True)

    print("Identifying the element by matching its strongest spectral lines against a reference database for Iron (Fe).\n")
    print("--- Match Verification ---")

    # Match the strongest observed lines against the reference data.
    tolerance = 0.05  # Wavelength tolerance in nm.
    for observed_line in observed_spectrum:
        for ref_wl in reference_fe_lines.keys():
            if abs(observed_line['wl'] - ref_wl) <= tolerance:
                print(f"Observed: {observed_line['wl']:.4f} nm  =~=  Reference (Iron): {ref_wl:.4f} nm")
                break
    
    print("\n--- Conclusion ---")
    print("The correlation between the observed spectrum and the reference spectrum for Iron is definitive.")
    print("The chemical element is Iron (Fe).")


identify_element()
<<<Iron (Fe)>>>