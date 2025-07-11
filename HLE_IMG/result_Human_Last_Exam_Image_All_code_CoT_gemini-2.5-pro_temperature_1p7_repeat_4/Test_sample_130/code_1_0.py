import sys

def identify_element_from_spectrum():
    """
    Identifies an element by comparing its observed spectral lines to a database of known spectra.
    """
    # Step 1: Database of known prominent spectral lines (in nanometers).
    SPECTRAL_DATA = {
        "Hydrogen": {"H-beta": 486.1, "H-alpha": 656.3},
        "Helium": {"blue-violet": 447.1, "cyan": 492.2, "green": 501.6, "yellow": 587.6, "red": 667.8},
        "Mercury": {"blue": 435.8, "green": 546.1, "yellow": 578.0}
    }

    # Step 2: Prominent lines estimated from the provided image.
    # The image shows strong lines in the blue, cyan, green, yellow, and red regions.
    observed_lines = {
        "blue": 447,
        "cyan": 492,
        "green": 502,
        "yellow": 588,
        "red": 668
    }
    
    print("Analyzing the spectral lines...")
    print("Observed prominent lines (estimated from image in nm):")
    for color, wavelength in observed_lines.items():
        print(f"- A bright {color} line at approximately {wavelength} nm")
    print("-" * 40)

    # Step 3: Comparison logic to find the best match.
    best_match_element = "Unknown"
    max_matches = 0
    best_match_details = {}
    tolerance = 5.0  # Allow a tolerance of +/- 5 nm for estimation.

    for element, known_lines in SPECTRAL_DATA.items():
        current_matches = 0
        match_details = {}
        for obs_color, obs_wl in observed_lines.items():
            for known_line_name, known_wl in known_lines.items():
                if abs(obs_wl - known_wl) <= tolerance:
                    current_matches += 1
                    match_details[f"Observed {obs_color} ({obs_wl} nm)"] = f"{element} ({known_wl} nm)"
                    break 
        
        if current_matches > max_matches:
            max_matches = current_matches
            best_match_element = element
            best_match_details = match_details

    # Step 4: Display the final result and the matching "equation".
    if best_match_element != "Unknown":
        print(f"Conclusion: The spectrum strongly matches the emission lines of {best_match_element}.")
        print("\nThe following comparisons confirm this identification:")
        
        # This part fulfills the "output each number in the final equation" requirement.
        for obs, known in best_match_details.items():
            obs_val = obs.split('(')[1].split(' ')[0]
            known_val = known.split('(')[1].split(' ')[0]
            print(f"Observed Line ({obs_val} nm) ≈ Known {best_match_element} Line ({known_val} nm)")
    else:
        print("Could not confidently identify the element from the provided data.")

identify_element_from_spectrum()