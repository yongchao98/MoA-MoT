import sys
import subprocess

# This script requires the 'colour-science' library.
# We will check if it's installed and provide instructions if it's not.
try:
    import numpy as np
    import colour
except ImportError:
    print("The 'colour-science' library is not installed.")
    print("Please install it by running the following command in your shell:")
    print("\npip install colour-science\n")
    sys.exit(1)

def check_romm_to_srgb_gamut():
    """
    Checks which ROMM RGB colors are outside the sRGB gamut.
    """
    # Conversion matrix from linear ProPhoto RGB (ROMM) to linear sRGB.
    # This matrix from the colour-science library accounts for the different
    # primaries and white points (D50 for ROMM, D65 for sRGB).
    M_romm_to_srgb = colour.RGB_COLOURSPACES['ProPhoto RGB'].matrix_ProPhotoRGB_to_sRGB

    # The colors to test, provided as non-linear ROMM RGB values.
    colors_to_test = {
        1: {'label': 'RGB(0, 0, 1)', 'value': np.array([0, 0, 1.0])},
        2: {'label': 'RGB(0, 1, 0)', 'value': np.array([0, 1.0, 0])},
        3: {'label': 'RGB(0, 0.5, 0.6)', 'value': np.array([0, 0.5, 0.6])},
        4: {'label': 'RGB(0.4, 0.5, 0.6)', 'value': np.array([0.4, 0.5, 0.6])},
        5: {'label': 'RGB(1, 1, 1)', 'value': np.array([1.0, 1.0, 1.0])}
    }

    out_of_gamut_indices = []

    print("Checking which ROMM RGB colors are out of the sRGB gamut...\n")

    for index, color_data in sorted(colors_to_test.items()):
        non_linear_romm = color_data['value']
        label = color_data['label']

        # Step 1: Linearize the ROMM RGB values.
        # ProPhoto RGB uses a gamma of 1.8 for its transfer function.
        linear_romm = colour.cctf_decoding(non_linear_romm, function='ProPhoto RGB')

        # Step 2: Convert from linear ROMM RGB to linear sRGB using the matrix.
        linear_srgb = np.dot(M_romm_to_srgb, linear_romm)

        # Step 3: Check if the resulting linear sRGB values are in the [0, 1] gamut.
        # A valid sRGB color must have all components >= 0 and <= 1.
        is_in_gamut = np.all((linear_srgb >= 0) & (linear_srgb <= 1))
        
        print(f"Checking color {index}) {label}:")
        # To make the output clear, we round the result for printing
        print(f"  Converted to linear sRGB: {np.round(linear_srgb, 4)}")

        if not is_in_gamut:
            print("  Result: Out of Gamut (cannot be represented by sRGB)\n")
            out_of_gamut_indices.append(index)
        else:
            print("  Result: In Gamut (can be represented by sRGB)\n")
            
    # Format the final answer
    if not out_of_gamut_indices:
        result_string = "none cannot"
    else:
        result_string = ", ".join(map(str, sorted(out_of_gamut_indices)))
    
    print("---")
    print("Final Answer:")
    print(f"The numbers of the colors that cannot be represented are: {result_string}")

if __name__ == "__main__":
    check_romm_to_srgb_gamut()