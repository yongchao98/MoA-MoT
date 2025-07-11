import matplotlib.colors as mcolors

def hex_to_rgb_linear(hex_color):
    """Converts a hex color string to a linear RGB tuple."""
    rgb_srgb = mcolors.to_rgb(hex_color)
    rgb_linear = []
    for val in rgb_srgb:
        if val <= 0.04045:
            rgb_linear.append(val / 12.92)
        else:
            rgb_linear.append(((val + 0.055) / 1.055) ** 2.4)
    return tuple(rgb_linear)

def calculate_luminance(hex_color):
    """Calculates the relative luminance of a hex color."""
    r, g, b = hex_to_rgb_linear(hex_color)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b

def analyze_palette(name, hex_codes):
    """Analyzes a palette and prints the luminance of each color."""
    print(f"Analyzing palette for {name}:")
    luminances = []
    for i, hex_code in enumerate(hex_codes):
        lum = calculate_luminance(hex_code)
        luminances.append(lum)
        print(f"  Color {i+1} ({hex_code}): Luminance = {lum:.4f}")
    
    # Check if luminances are distinct
    # Round to 2 decimal places for practical distinctness check
    rounded_lums = [round(l, 2) for l in luminances]
    if len(set(rounded_lums)) == len(rounded_lums):
        print("  Result: All colors have distinct luminance values. GOOD for monochromatic vision.\n")
    else:
        print("  Result: Some colors have similar or identical luminance. BAD for monochromatic vision.\n")

# --- Palette Definitions ---

# Palette from Plot 2: pals::ocean.balance(5)
# This is a diverging palette designed for perceptual uniformity.
pal1_ocean_balance = ["#4B3A71", "#8B84A3", "#D5D5D5", "#9DB8C4", "#2C7396"]

# Palette from Plot 6: scales::hue_pal()(5) (ggplot2 default)
# This palette has colors of similar luminance.
pal5_ggplot_default = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

# --- Analysis ---

print("Monochromatic vision relies on luminance (brightness) to distinguish colors.")
print("A suitable palette must contain colors with clearly different luminance values.\n")

analyze_palette("Plot 2", pal1_ocean_balance)
analyze_palette("Plot 6 (ggplot2 default)", pal5_ggplot_default)

print("Analysis of other plots based on their generation method:")
print("- Plots 3 & 4: Generated with a constant lightness parameter (l=60). All colors will have the same luminance.")
print("- Plot 5: Generated with lightness values sampled from (20, 50, 70, 20, 90). The repeated '20' means two colors will have identical luminance.")

print("\nConclusion: Only Plot 2 has a palette where all colors are clearly distinguishable by luminance.")
print("Therefore, the interpretable plot for someone with full monochromatic vision is Plot 2.")
