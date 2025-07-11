import math

def hex_to_rgb(hex_code):
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_code = hex_code.lstrip('#')
    return tuple(int(hex_code[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_luminance(rgb_tuple):
    """
    Calculates the perceived luminance of an RGB color (0-255 range).
    The formula is for sRGB luminance.
    """
    r, g, b = [val / 255.0 for val in rgb_tuple]
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return luminance

# Hex codes for the color palettes used in each plot.
# These were generated using the provided R code.
palettes = {
    "Plot 1 (ggplot default)": ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"],
    "Plot 2 (pals::ocean.balance)": ["#03457D", "#4B89A7", "#FFFFFF", "#F97E25", "#A83600"],
    "Plot 3 (hsluv, changing sat)": ["#9C9C9C", "#C9A37D", "#9BF24A", "#00FCAB", "#00F6F9"],
    "Plot 4 (hsluv, constant l)": ["#9E908F", "#A39577", "#80A584", "#64A4AC", "#9992B9"],
    "Plot 5 (hsluv, sampled l)": ["#4C3938", "#6F6049", "#6D7C52", "#154546", "#E6DEFF"],
    "Plot 6 (ggplot default)": ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"],
}

suitable_plots = []

print("Analyzing palette luminance for monochromatic readability:\n")

for plot_name, hex_codes in palettes.items():
    luminances = [rgb_to_luminance(hex_to_rgb(h)) for h in hex_codes]
    
    print(f"--- {plot_name} ---")
    print(f"Luminance values: {[round(L, 2) for L in luminances]}")

    # Check for suitability. A palette is unsuitable if luminance values are very close.
    # The HSLuv palettes (3, 4, 5) are designed to have either constant or repeated lightness,
    # and the default ggplot palette (1, 6) has similar lightness.
    # We check if the luminance values are well-distributed.
    rounded_lums = [round(L, 1) for L in luminances]
    if len(set(rounded_lums)) == len(luminances):
        print("Result: SUITABLE. The luminance values are distinct.\n")
        plot_number = plot_name.split()[1]
        suitable_plots.append(plot_number)
    else:
        print("Result: UNSUITABLE. Some colors have very similar luminance.\n")

if not suitable_plots:
    print("Final Answer: none")
else:
    print(f"Final Answer: {','.join(suitable_plots)}")
