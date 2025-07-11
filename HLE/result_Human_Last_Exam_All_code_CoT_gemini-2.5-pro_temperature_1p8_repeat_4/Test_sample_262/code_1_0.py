def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    """Converts a hex color string to an (R, G, B) tuple (0-255)."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def calculate_luminance(rgb: tuple[int, int, int]) -> float:
    """Calculates the relative luminance (0-100) for an (R, G, B) tuple."""
    srgb = [val / 255.0 for val in rgb]
    # Linearize sRGB values
    rgb_linear = [
        (c / 12.92) if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
        for c in srgb
    ]
    # Calculate luminance using the standard formula and scale to 0-100
    return (0.2126 * rgb_linear[0] + 0.7152 * rgb_linear[1] + 0.0722 * rgb_linear[2]) * 100

def is_palette_monochromatic_friendly(luminances: list[float], threshold: float = 4.0) -> bool:
    """Checks if luminance values are distinct enough for monochromatic vision."""
    if len(luminances) < 2:
        return True
    sorted_lum = sorted(luminances)
    for i in range(len(sorted_lum) - 1):
        if (sorted_lum[i+1] - sorted_lum[i]) < threshold:
            return False
    return True

# --- Analysis ---

# Palettes used in the plots
# Plot 1 & 6: Default ggplot2 palette (scales::hue_pal). Varies hue, not ideal for luminance contrast.
palette_plot1_6 = ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']

# Plot 2: A diverging palette (pals::ocean.balance). Designed to vary lightness.
palette_plot2 = ['#0D47A1', '#42A5F5', '#FFFFFF', '#EF5350', '#B71C1C']

# Plots 3, 4, 5 are unsuitable by their programmatic definition in R:
# Plot 3 (pal2): Uses HSLuv color space with varying Saturation but fixed Lightness (L=60).
#               All colors will have the same luminance and be indistinguishable in grayscale.
# Plot 4 (pal3): Uses HSLuv color space with fixed Saturation and Lightness (L=60).
#               Same issue as Plot 3.
# Plot 5 (pal4): Uses HSLuv color space with Lightness sampled from a set containing duplicates {20,50,70,20,90}.
#               Two of the five categories will have identical luminance (20), making them indistinguishable.

# We only need to numerically analyze the palettes for Plots 1, 2, and 6.
palettes_to_check = {
    1: palette_plot1_6,
    2: palette_plot2,
    6: palette_plot1_6,
}

suitable_plots = []
print("Analyzing palettes for monochromatic interpretability...\n")

for plot_num, palette_hex in palettes_to_check.items():
    print(f"--- Analyzing Plot {plot_num} ---")
    luminances = [calculate_luminance(hex_to_rgb(c)) for c in palette_hex]
    is_friendly = is_palette_monochromatic_friendly(luminances)

    print(f"Luminance values: {[round(l, 1) for l in sorted(luminances)]}")
    if is_friendly:
        print("Result: SUITABLE. The luminance values are clearly distinct.")
        if plot_num not in suitable_plots:
             suitable_plots.append(plot_num)
    else:
        print("Result: NOT SUITABLE. Some colors have very similar luminance levels.")
    print("-" * 25)

print("\n--- Conclusion ---")
print("Plots 3, 4, and 5 are unsuitable by design, as they produce colors with identical or duplicate luminance values.")
if not suitable_plots:
    final_answer = "none"
    print("None of the tested plots are suitable.")
else:
    final_answer = ",".join(map(str, sorted(suitable_plots)))
    print(f"The only plot that uses a color palette suitable for monochromatic vision is: Plot {final_answer}")

print(f"\n<<<{final_answer}>>>")