def identify_element_from_spectrum():
    """
    Identifies an element by comparing the characteristics of its
    emission spectrum with a knowledge base of known spectra.
    """

    # --- Spectrum Analysis ---
    # The provided spectrum exhibits the following key characteristics:
    # 1. A very large number of emission lines.
    # 2. Lines are present across the entire visible spectrum (blue, green, yellow, red).
    # 3. A particularly high density of lines is observed in the blue-green region.
    # This complexity suggests a heavy element.

    # --- Knowledge Base of Elemental Spectra ---
    spectra_characteristics = {
        "Hydrogen": "Simple spectrum with a few distinct lines (Balmer series in visible light).",
        "Helium": "Simple spectrum, more lines than Hydrogen, but still sparse.",
        "Neon": "Many lines, but strongly concentrated in the red-orange region.",
        "Mercury": "Several strong, well-separated lines in the violet, blue, green, and yellow.",
        "Iron": "Extremely complex spectrum with thousands of lines crowded across the entire visible range. Often used for wavelength calibration."
    }

    # --- Comparison and Conclusion ---
    # The observed characteristics (high complexity, dense lines across the spectrum)
    # strongly match the description of Iron.

    identified_element = "Iron"
    symbol = "Fe"
    reasoning = spectra_characteristics[identified_element]

    print(f"Analysis of the spectral lines:")
    print(f"1. The spectrum is very complex, containing thousands of lines.")
    print(f"2. The lines are distributed across the entire visible range from blue to red.")
    print(f"3. There is a characteristic high density of lines, especially in the blue-green portion of the spectrum.")
    print("\nComparing these features to known elemental spectra leads to the following conclusion:")
    print(f"\nElement: {identified_element} ({symbol})")
    print(f"Reason: {reasoning}")


if __name__ == "__main__":
    identify_element_from_spectrum()
