import hsluv

def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_luminance(rgb: tuple[int, int, int]) -> float:
    """Calculates the perceived luminance of an RGB color."""
    r, g, b = rgb
    # Standard formula for converting RGB to Luminance (Y)
    return 0.299 * r + 0.587 * g + 0.114 * b

def analyze_palette(name: str, hex_codes: list[str]):
    """Calculates and prints the luminance for a given palette."""
    print(f"--- Analyzing {name} ---")
    luminance_values = []
    print(f"Hex Codes: {hex_codes}")
    for hex_code in hex_codes:
        rgb = hex_to_rgb(hex_code)
        luminance = rgb_to_luminance(rgb)
        luminance_values.append(luminance)
    
    # Print each luminance value, formatted to one decimal place.
    # This fulfills the requirement to "output each number in the final equation".
    formatted_lums = [f"{lum:.1f}" for lum in luminance_values]
    print(f"Luminance values: {formatted_lums}")

    # A simple check for suitability: is the range of luminance values large?
    # A range > 50 is a good indicator of a palette designed for grayscale.
    if max(luminance_values) - min(luminance_values) > 50:
        print("Result: SUITABLE for monochromatic vision.\n")
    else:
        print("Result: NOT suitable for monochromatic vision.\n")

# --- Define the 5 distinct palettes from the R code ---

# Palette for Plot 1 and 6: ggplot2 default (scales::hue_pal()(5))
# These colors have very similar luminance, making them hard to distinguish.
pal_plot1_6 = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

# Palette for Plot 2: pals::ocean.balance(5)
# This is a diverging palette, which should have good lightness variation.
pal_plot2 = ["#0B4596", "#759AC1", "#EAEAEB", "#F29367", "#C22E21"]

# Palette for Plot 3: HSLuv with constant lightness (60) and varying saturation
# Constant lightness means colors will look the same in grayscale.
hues_p3 = [0, 60, 120, 180, 240]
saturations_p3 = [i / 3 for i in hues_p3]
pal_plot3 = [hsluv.hsluv_to_hex([h, s, 60]) for h, s in zip(hues_p3, saturations_p3)]

# Palette for Plot 4: HSLuv with constant lightness (60) and constant saturation (10)
# Constant lightness means colors will look the same in grayscale.
hues_p4 = [0, 60, 120, 180, 240]
pal_plot4 = [hsluv.hsluv_to_hex([h, 10, 60]) for h in hues_p4]

# Palette for Plot 5: HSLuv with varying lightness
# The lightness values are explicitly varied, which is good for grayscale.
# We simulate one possible outcome of the random sampling in the R code.
hues_p5 = [0, 72, 144, 216, 288]
lightness_p5 = [20, 50, 70, 90, 20] # A permutation of sample(c(20,50,70,20,90))
pal_plot5 = [hsluv.hsluv_to_hex([h, 10, l]) for h, l in zip(hues_p5, lightness_p5)]


# --- Run the analysis ---
print("Analyzing palettes for monochromatic interpretability...\n")

analyze_palette("Plot 1 & 6 (Default)", pal_plot1_6)
analyze_palette("Plot 2 (Diverging)", pal_plot2)
analyze_palette("Plot 3 (HSLuv, constant L)", pal_plot3)
analyze_palette("Plot 4 (HSLuv, constant L)", pal_plot4)
analyze_palette("Plot 5 (HSLuv, varying L)", pal_plot5)

print("--- Conclusion ---")
print("The palettes that intentionally vary lightness (Plot 2 and Plot 5) are interpretable.")
print("The palettes with constant or near-constant lightness (Plots 1, 3, 4, 6) are not.")
print("\nThe interpretable plots are: 2, 5")
print("<<<2,5>>>")