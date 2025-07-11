# Required libraries: colour-science and numpy
# You may need to install them first, for example: pip install colour-science numpy
import numpy as np
import colour

# --- Step 1: Define the list of ROMM RGB colors to be tested ---
colors_to_test = {
    1: np.array([0, 0, 1]),
    2: np.array([0, 1, 0]),
    3: np.array([0, 0.5, 0.6]),
    4: np.array([0.4, 0.5, 0.6]),
    5: np.array([1, 1, 1])
}

# Define the input and output color spaces for the conversion function
# ROMM RGB is officially known as 'ProPhoto RGB'
INPUT_COLOURSPACE = 'ProPhoto RGB'
OUTPUT_COLOURSPACE = 'sRGB'

out_of_gamut_indices = []

print(f"Analyzing which ROMM RGB colors are outside the {OUTPUT_COLOURSPACE} gamut...\n")

# --- Step 2 & 3: Iterate through colors, convert them, and check their gamut ---
for index, romm_rgb in colors_to_test.items():
    
    # Convert the color from ROMM RGB to sRGB
    srgb_color = colour.RGB_convert(romm_rgb, INPUT_COLOURSPACE, OUTPUT_COLOURSPACE)
    
    # Check if any component is outside the valid [0.0, 1.0] range
    is_out_of_gamut = np.any((srgb_color < 0.0) | (srgb_color > 1.0))
    
    status = "Out of gamut (Cannot be represented)" if is_out_of_gamut else "In gamut (Can be represented)"
    # Round the sRGB values for cleaner display
    srgb_rounded = np.round(srgb_color, 4)
    print(f"Color #{index}: ROMM RGB {romm_rgb} -> sRGB {srgb_rounded} -> {status}")
    
    # --- Step 4: Collect the numbers of out-of-gamut colors ---
    if is_out_of_gamut:
        out_of_gamut_indices.append(index)

# --- Step 5: Sort the final list and print the result ---
out_of_gamut_indices.sort()

if not out_of_gamut_indices:
    final_answer_string = "none cannot"
else:
    final_answer_string = ", ".join(map(str, out_of_gamut_indices))

print("\n---------------------------------------------------------------------")
print("Final Answer: The numbers of the colors that cannot be represented are:")
print(final_answer_string)
print("---------------------------------------------------------------------")
