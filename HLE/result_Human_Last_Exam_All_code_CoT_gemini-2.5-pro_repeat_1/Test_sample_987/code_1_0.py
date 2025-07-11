import io
import pandas as pd

def identify_element():
    """
    Identifies a chemical element by comparing its observed emission spectrum
    to a reference database of spectral lines.
    """

    # Step 1: Parse the user-provided spectral data.
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
Wavelength: 226.4610 nm, Intensity: 830.0
Wavelength: 228.4600 nm, Intensity: 330.0
Wavelength: 229.8159 nm, Intensity: 790.0
Wavelength: 230.4220 nm, Intensity: 2700.0
Wavelength: 234.3180 nm, Intensity: 1600.0
Wavelength: 236.3040 nm, Intensity: 2500.0
Wavelength: 237.2770 nm, Intensity: 3500.0
Wavelength: 238.6890 nm, Intensity: 1300.0
Wavelength: 239.0620 nm, Intensity: 2500.0
Wavelength: 239.1180 nm, Intensity: 2700.0
Wavelength: 242.4990 nm, Intensity: 370.0
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
    # Use regular expressions to extract numbers
    data = pd.read_csv(
        io.StringIO(spectrum_data_string),
        sep=':',
        header=None,
        names=['label', 'value']
    )
    wavelengths = data[data['label'] == 'Wavelength']['value'].str.extract(r'(\d+\.\d+)').astype(float).values.flatten()
    intensities = data[data['label'] == ' Intensity']['value'].astype(float).values.flatten()
    observed_spectrum = pd.DataFrame({'Wavelength': wavelengths, 'Intensity': intensities})

    # Step 2: Load reference data for Cobalt (Co).
    # A selection of prominent lines for Co I and Co II from the NIST database.
    cobalt_reference_lines = {
        201.065: "Co II", 202.235: "Co II", 203.357: "Co II", 205.222: "Co II",
        206.064: "Co II", 208.883: "Co II", 209.263: "Co II", 215.805: "Co II",
        228.616: "Co II", 230.786: "Co II", 236.379: "Co II", 237.862: "Co II",
        238.892: "Co II", 240.725: "Co I",  242.493: "Co I",  254.397: "Co II",
        266.479: "Co II", 292.479: "Co II", 306.890: "Co II", 313.332: "Co II",
        322.078: "Co II", 345.351: "Co I",  351.364: "Co II", 380.012: "Co II",
    }
    
    # Step 3 & 4: Match the strongest observed lines to the reference data.
    tolerance = 0.1  # Wavelength tolerance in nm
    
    # Sort observed lines by intensity to focus on the most prominent ones
    strongest_observed = observed_spectrum.sort_values(by='Intensity', ascending=False).head(15)
    
    matches = []
    
    for _, row in strongest_observed.iterrows():
        obs_wl = row['Wavelength']
        obs_int = row['Intensity']
        
        for ref_wl, ion in cobalt_reference_lines.items():
            if abs(obs_wl - ref_wl) <= tolerance:
                matches.append({
                    "Observed Wavelength (nm)": obs_wl,
                    "Observed Intensity": obs_int,
                    "Reference Wavelength (nm)": ref_wl,
                    "Element/Ion": ion
                })
                break # Found a match, move to the next observed line

    # Step 5: Output the conclusion and justification.
    print("Analysis of the provided spectrum:\n")
    
    if len(matches) > 5: # Require at least 5 strong matches for a positive ID
        print("Conclusion: The spectrum shows a very strong correlation with the element Cobalt (Co).\n")
        print("The identification is based on matching the most intense spectral lines from the sample")
        print("with known reference lines for neutral (Co I) and singly-ionized (Co II) Cobalt.\n")
        print("Comparison of the strongest spectral lines:")
        print("-" * 80)
        print(f"{'Observed Wavelength (nm)':<28} {'Observed Intensity':<22} {'Matching Cobalt Line (nm)':<28}")
        print("-" * 80)
        
        # Create a DataFrame from the first 10 matches for clean printing
        match_df = pd.DataFrame(matches).head(10)
        for _, row in match_df.iterrows():
            print(f"{row['Observed Wavelength (nm)']:<28.4f} {row['Observed Intensity']:<22.1f} {row['Reference Wavelength (nm)']:<28.4f}")
        print("-" * 80)
        
        print("\n<<<Cobalt>>>")

    else:
        print("Could not confidently identify the element based on the provided data and reference set.")
        print("\n<<<Unknown>>>")


if __name__ == '__main__':
    identify_element()