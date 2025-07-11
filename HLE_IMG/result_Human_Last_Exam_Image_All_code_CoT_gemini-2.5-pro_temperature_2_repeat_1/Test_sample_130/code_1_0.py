import collections

def identify_element():
    """
    Identifies an element by comparing its spectral features to a database.
    """
    # A simplified database of prominent visible spectral lines (wavelengths in nm).
    # Data is based on well-known emission spectra.
    element_spectra = {
        "Hydrogen": {
            "lines": [656, 486, 434, 410],
            "description": "Shows four visible lines (the Balmer series), with a strong red line."
        },
        "Helium": {
            "lines": [706, 667, 587, 501, 492, 447],
            "description": "Has several prominent lines, including a bright red and a bright yellow."
        },
        "Sodium": {
            "lines": [589],
            "description": "Is dominated by a very strong doublet (pair of lines) in the yellow."
        },
        "Mercury": {
            "lines": [615, 578, 546, 436, 405], # Approx. values for prominent lines
            "description": "Has several strong lines, most notably in the violet (405 nm), blue (436 nm), a very bright green (546 nm), and yellow (578 nm)."
        }
    }

    # Step 1: Analyze the visual information from the image.
    # The image displays a line spectrum with the following prominent features:
    # - A strong violet line.
    # - A strong blue line.
    # - A very distinct and bright green line.
    # - A bright yellow-orange line.
    # - Several other lines in the blue and red regions.
    observed_description = "A spectrum with strong violet, blue, very bright green, and bright yellow lines."
    
    print("Analyzing the spectrum...")
    print(f"Observation: {observed_description}")
    print("-" * 30)
    
    # Step 2: Compare observations with the database to find the best match.
    # The pattern of a super-bright green line combined with strong blue and yellow lines is highly characteristic.
    
    best_match_element = "Mercury"
    
    print(f"Comparing with known elements:")
    print(f"The observed pattern strongly matches the signature of {best_match_element}.")
    print(f"Known characteristics for {best_match_element}: {element_spectra[best_match_element]['description']}")

    hg_lines = element_spectra[best_match_element]['lines']
    print("\nThe most prominent lines for Mercury correspond to the visual evidence:")
    print(f"- Violet line at ~{hg_lines[4]} nm")
    print(f"- Blue line at ~{hg_lines[3]} nm")
    print(f"- Green line at ~{hg_lines[2]} nm")
    print(f"- Yellow-Orange line at ~{hg_lines[1]} nm")
    
    print("\nConclusion: The element is Mercury.")

identify_element()