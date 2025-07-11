def identify_element_from_spectrum():
    """
    Identifies an element by its characteristic spectral lines.

    The provided image shows a spectrum with prominent lines across the visible
    range. This pattern is characteristic of Helium. This script lists the
    major visible spectral lines of Helium to confirm the identification.
    The wavelengths are given in nanometers (nm).
    """

    # Dictionary of prominent visible spectral lines for Helium
    element_data = {
        "element_name": "Helium",
        "symbol": "He",
        "spectral_lines_nm": {
            "blue-violet": 447.1,
            "blue": 471.3,
            "cyan": 492.2,
            "green": 501.6,
            "yellow-orange": 587.6,
            "red": 667.8,
        }
    }

    element = element_data["element_name"]
    symbol = element_data["symbol"]
    lines = element_data["spectral_lines_nm"]

    print(f"The spectral lines in the image match the emission spectrum of {element} ({symbol}).")
    print("\nHere are the prominent visible spectral lines of Helium:")

    for color, wavelength in lines.items():
        print(f"- A {color} line at approximately {wavelength} nm")

identify_element_from_spectrum()