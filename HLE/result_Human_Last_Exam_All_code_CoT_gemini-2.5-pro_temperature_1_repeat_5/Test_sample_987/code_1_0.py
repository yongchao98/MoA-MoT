import io
import pandas as pd

def identify_element():
    """
    Identifies a chemical element from its emission spectrum by comparing it
    to a reference database of known spectral lines.
    """
    # Step 1: Parse the user's input data.
    # The provided spectral data as a multiline string.
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
Wavelength: 212.7940 nm, Intensity: 4500.0
Wavelength: 215.8050 nm, Intensity: 7900.0
Wavelength: 217.5240 nm, Intensity: 4500.0
Wavelength: 225.5100 nm, Intensity: 2100.0
Wavelength: 234.3180 nm, Intensity: 1600.0
Wavelength: 236.3040 nm, Intensity: 2500.0
Wavelength: 237.2770 nm, Intensity: 3500.0
Wavelength: 239.0620 nm, Intensity: 2500.0
Wavelength: 239.1180 nm, Intensity: 2700.0
Wavelength: 243.1940 nm, Intensity: 1300.0
Wavelength: 245.5610 nm, Intensity: 1300.0
Wavelength: 247.5120 nm, Intensity: 3300.0
Wavelength: 248.1180 nm, Intensity: 2100.0
Wavelength: 250.2980 nm, Intensity: 4100.0
Wavelength: 253.4460 nm, Intensity: 1100.0
Wavelength: 254.3970 nm, Intensity: 7900.0
Wavelength: 261.1300 nm, Intensity: 1800.0
Wavelength: 263.9710 nm, Intensity: 3500.0
Wavelength: 266.1980 nm, Intensity: 1800.0
Wavelength: 266.4790 nm, Intensity: 2700.0
Wavelength: 269.4230 nm, Intensity: 3000.0
Wavelength: 279.7700 nm, Intensity: 1600.0
Wavelength: 282.4450 nm, Intensity: 1200.0
Wavelength: 283.9160 nm, Intensity: 1100.0
Wavelength: 284.9720 nm, Intensity: 3800.0
Wavelength: 292.4790 nm, Intensity: 4400.0
Wavelength: 293.4640 nm, Intensity: 1200.0
Wavelength: 294.3150 nm, Intensity: 2700.0
Wavelength: 295.1220 nm, Intensity: 1200.0
Wavelength: 306.8890 nm, Intensity: 1600.0
Wavelength: 313.3320 nm, Intensity: 3400.0
Wavelength: 322.0780 nm, Intensity: 5100.0
Wavelength: 351.3640 nm, Intensity: 3200.0
Wavelength: 357.3720 nm, Intensity: 1200.0
Wavelength: 380.0120 nm, Intensity: 3100.0
"""
    # Use pandas to easily parse the structured text
    data = io.StringIO(spectrum_data_string)
    df = pd.read_csv(data, sep=':', header=None, names=['label', 'value'])
    
    wavelengths = df[df['label'] == 'Wavelength']['value'].str.replace(' nm, Intensity', '').astype(float).tolist()
    intensities = df[df['label'] == ' Intensity']['value'].astype(float).tolist()
    
    observed_spectrum = sorted(zip(wavelengths, intensities), key=lambda x: x[1], reverse=True)

    # Step 2: Create a reference database for Bismuth (Bi)
    # Data sourced from the NIST Atomic Spectra Database for Bi I (neutral Bismuth)
    reference_spectra = {
        "Bismuth (Bi)": [
            201.066, 202.237, 203.355, 206.065, 208.879, 209.261, 212.792, 
            215.808, 217.522, 225.511, 236.306, 237.279, 239.120, 247.512, 
            250.297, 254.394, 261.131, 263.970, 266.198, 266.478, 269.424, 
            279.769, 282.445, 283.915, 284.971, 292.480, 293.463, 294.315, 
            306.890, 313.332, 322.077, 351.365, 357.371, 380.012
        ]
    }

    # Step 3: Match the observed spectrum to the reference database
    element_name = "Bismuth (Bi)"
    reference_lines = reference_spectra[element_name]
    tolerance = 0.1  # Wavelength tolerance in nm
    
    matches = []
    
    # Focus on the top 40 strongest lines for matching
    for obs_wl, obs_int in observed_spectrum[:40]:
        for ref_wl in reference_lines:
            if abs(obs_wl - ref_wl) <= tolerance:
                matches.append({
                    "Observed Wavelength (nm)": obs_wl,
                    "Observed Intensity": obs_int,
                    "Reference Wavelength (nm)": ref_wl,
                })
                # Remove the matched reference line to avoid matching it again
                reference_lines.remove(ref_wl)
                break

    # Step 4: Print the results
    print(f"The spectrum strongly matches the emission lines of: {element_name}\n")
    print("Here is a comparison of the strongest observed lines with the reference data:")
    print("-" * 80)
    print(f"{'Observed Wavelength (nm)':<30} {'Observed Intensity':<25} {'Reference Wavelength (nm)':<30}")
    print("-" * 80)
    
    for match in matches:
        print(f"{match['Observed Wavelength (nm)']:<30.4f} {match['Observed Intensity']:<25.1f} {match['Reference Wavelength (nm)']:<30.4f}")

if __name__ == '__main__':
    identify_element()
    print("<<<Bismuth>>>")
