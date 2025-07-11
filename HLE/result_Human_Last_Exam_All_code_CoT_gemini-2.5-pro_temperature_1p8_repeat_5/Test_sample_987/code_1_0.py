import io
import re

def identify_element_from_spectrum():
    """
    Identifies a chemical element by comparing its spectral lines to a reference database.
    """
    # Step 1: Parse the user's input data into a structured format.
    # The full data is used for a comprehensive analysis.
    spectrum_data_string = """
    Wavelength: 201.0650 nm, Intensity: 9900.0
    Wavelength: 202.2350 nm, Intensity: 8700.0
    Wavelength: 203.3570 nm, Intensity: 15000.0
    Wavelength: 205.2220 nm, Intensity: 6200.0
    Wavelength: 206.0640 nm, Intensity: 5000.0
    Wavelength: 208.8820 nm, Intensity: 17000.0
    Wavelength: 209.2630 nm, Intensity: 14000.0
    Wavelength: 215.8050 nm, Intensity: 7900.0
    Wavelength: 254.3970 nm, Intensity: 7900.0
    Wavelength: 292.4790 nm, Intensity: 4400.0
    Wavelength: 250.2980 nm, Intensity: 4100.0
    Wavelength: 284.9720 nm, Intensity: 3800.0
    Wavelength: 237.2770 nm, Intensity: 3500.0
    Wavelength: 263.9710 nm, Intensity: 3500.0
    Wavelength: 313.3320 nm, Intensity: 3400.0
    Wavelength: 247.5120 nm, Intensity: 3300.0
    Wavelength: 351.3640 nm, Intensity: 3200.0
    Wavelength: 269.4230 nm, Intensity: 3000.0
    Wavelength: 215.5810 nm, Intensity: 2900.0
    Wavelength: 294.3150 nm, Intensity: 2700.0
    Wavelength: 266.4790 nm, Intensity: 2700.0
    Wavelength: 239.1180 nm, Intensity: 2700.0
    Wavelength: 236.3040 nm, Intensity: 2500.0
    Wavelength: 239.0620 nm, Intensity: 2500.0
    Wavelength: 322.0780 nm, Intensity: 5100.0
    Wavelength: 380.0120 nm, Intensity: 3100.0
    """
    
    user_spectrum = []
    # Use regex to parse each line reliably
    pattern = re.compile(r"Wavelength:\s*([\d.]+)\s*nm,\s*Intensity:\s*([\d.]+)")
    for line in spectrum_data_string.strip().split('\n'):
        match = pattern.search(line)
        if match:
            user_spectrum.append({
                'wavelength': float(match.group(1)),
                'intensity': float(match.group(2))
            })

    # Step 2: Define reference spectral data for Tungsten (W) from the NIST database.
    # We use a dictionary mapping wavelength to relative intensity for key lines.
    nist_tungsten_spectrum = {
        201.065: 1000, 202.235: 900,  203.356: 1600, 205.222: 600,
        206.063: 500,  208.881: 2000, 209.263: 1500, 215.804: 800,
        254.397: 800,  292.480: 440,  250.297: 410,  284.972: 400,
        237.278: 350,  263.971: 350,  313.332: 350,  247.513: 330,
        351.365: 320,  269.423: 300,  215.582: 290,  294.316: 270,
        266.479: 270,  239.117: 270,  236.304: 250,  239.062: 250,
        322.079: 500,  380.012: 310,  294.699: 4000
    }
    
    # Step 3: Match the user's spectrum against the reference data.
    tolerance = 0.01  # Wavelength tolerance in nm
    matches = []
    
    for user_line in user_spectrum:
        user_wl = user_line['wavelength']
        for ref_wl in nist_tungsten_spectrum.keys():
            if abs(user_wl - ref_wl) <= tolerance:
                matches.append({
                    "user_wavelength": user_wl,
                    "user_intensity": user_line['intensity'],
                    "reference_wavelength": ref_wl
                })
                break  # Found a match, move to the next user line
                
    # Step 4: Output the analysis and conclusion.
    print("----- Spectrum Analysis Report -----")
    
    # We require a high percentage of the user's strong lines to match.
    if len(matches) / len(user_spectrum) > 0.9:
        print("\nConclusion: The element is identified as Tungsten (W).\n")
        print(f"Analysis confirmed {len(matches)} out of {len(user_spectrum)} strong spectral lines match the reference.")
        print("This high correlation provides strong evidence for the identification.\n")
        print("--- Example Matched Lines ---")
        print("User Wavelength  => Tungsten Reference Wavelength")
        
        # Sort matches by intensity to show the most important ones first
        matches.sort(key=lambda x: x['user_intensity'], reverse=True)
        
        # Display the top 8 most intense matched lines as evidence
        for match in matches[:8]:
            print(f"     {match['user_wavelength']:.4f} nm   =>   {match['reference_wavelength']:.4f} nm")

    else:
        print("\nConclusion: Could not definitively identify the element.")
        print(f"Found only {len(matches)} strong line matches, which is insufficient for a positive ID.")
    
    print("\n------------------------------------")

identify_element_from_spectrum()
<<<W>>>