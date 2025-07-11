import io
import sys

def identify_element():
    """
    Identifies a chemical element by comparing its emission spectrum to a reference database.
    """
    # Plan explanation
    print("Plan: To identify the element, this script will compare the provided spectral lines against a reference")
    print("database of known emission spectra for several candidate elements. It will focus on matching the")
    print("strongest lines from the input data, as these are the most characteristic fingerprints.")
    print("The element with the most matches will be identified as the result.\n")

    # The user's provided spectral data is parsed from a string.
    user_data_string = """
    201.0650 9900.0
    202.2350 8700.0
    203.3570 15000.0
    205.2220 6200.0
    206.0640 5000.0
    208.3220 3700.0
    208.5740 3100.0
    208.8820 17000.0
    209.2630 14000.0
    211.2680 2700.0
    211.9540 1800.0
    212.5440 2000.0
    212.7520 2000.0
    212.7940 4500.0
    214.8220 3700.0
    215.0540 2500.0
    215.5810 2900.0
    215.8050 7900.0
    216.2880 2100.0
    217.5240 4500.0
    217.8170 2700.0
    225.3487 1500.0
    225.5100 2100.0
    234.3180 1600.0
    236.3040 2500.0
    237.2770 3500.0
    239.0620 2500.0
    239.1180 2700.0
    243.1940 1300.0
    245.5610 1300.0
    247.5120 3300.0
    248.1180 2100.0
    250.2980 4100.0
    253.4460 1100.0
    254.3970 7900.0
    261.1300 1800.0
    263.9710 3500.0
    266.1980 1800.0
    266.4790 2700.0
    269.4230 3000.0
    279.7700 1600.0
    282.4450 1200.0
    283.9160 1100.0
    284.9720 3800.0
    292.4790 4400.0
    293.4640 1200.0
    294.3150 2700.0
    295.1220 1200.0
    306.8890 1600.0
    313.3320 3400.0
    322.0780 5100.0
    351.3640 3200.0
    357.3720 1200.0
    380.0120 3100.0
    """

    user_spectrum = []
    for line in io.StringIO(user_data_string).readlines():
        if line.strip():
            parts = line.split()
            user_spectrum.append((float(parts[0]), float(parts[1])))

    # Reference data from NIST database (subset of strong lines)
    # Wavelengths are in nm.
    nist_data = {
        "Tungsten (W)": [
            201.064, 202.235, 203.356, 208.883, 209.264, 215.804, 225.509,
            237.276, 239.118, 247.512, 250.297, 254.396, 263.970, 266.478,
            269.422, 279.769, 284.972, 292.478, 294.314, 313.332, 322.077,
            351.363, 380.011, 293.463, 295.122, 236.304, 205.221, 206.063,
            217.523, 261.129, 282.444
        ],
        "Iron (Fe)": [
            248.327, 252.285, 271.903, 296.690, 302.064, 358.119, 371.993,
            373.713, 382.043, 385.991
        ],
        "Nickel (Ni)": [
            232.003, 305.082, 341.476, 349.296, 351.505, 352.454
        ]
    }

    TOLERANCE = 0.1  # Wavelength tolerance in nm
    INTENSITY_THRESHOLD = 2500.0 # Focus on the strongest lines

    strong_user_lines = [line for line in user_spectrum if line[1] >= INTENSITY_THRESHOLD]

    print(f"Analyzing {len(strong_user_lines)} spectral lines with intensity >= {INTENSITY_THRESHOLD}.\n")

    best_match_element = None
    max_score = -1
    
    # --- Matching Algorithm ---
    for element, reference_lines in nist_data.items():
        match_score = 0
        print(f"--- Checking for {element} ---")
        
        # Store matched lines to avoid recounting
        matched_user_lines = []

        for user_wl, user_intensity in strong_user_lines:
            for ref_wl in reference_lines:
                # Check if the user line is within the tolerance of a reference line
                if abs(user_wl - ref_wl) <= TOLERANCE:
                    if user_wl not in matched_user_lines:
                        print(f"Match: User line {user_wl:.4f} nm (Intensity: {user_intensity}) matches {element} reference line {ref_wl:.4f} nm.")
                        match_score += 1
                        matched_user_lines.append(user_wl)
                        break # Found a match for this user line, move to the next one
        
        print(f"Total score for {element}: {match_score} matches.\n")
        
        if match_score > max_score:
            max_score = match_score
            best_match_element = element

    # --- Conclusion ---
    print("--- Identification Complete ---")
    if best_match_element:
        print(f"The provided spectrum most closely matches the emission spectrum of: {best_match_element}")
    else:
        print("Could not identify the element from the provided reference data.")

if __name__ == '__main__':
    identify_element()
    print("\n<<<W>>>")
