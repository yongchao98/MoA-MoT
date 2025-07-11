import math

def hex_to_rgb_normalized(hex_color):
    """Converts a hex color string (e.g., '#RRGGBB') to normalized (r, g, b) values."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))

def calculate_luminance(rgb_normalized):
    """Calculates luminance on a 0-100 scale from normalized RGB."""
    r, g, b = rgb_normalized
    # Standard formula for relative luminance (Y)
    return (0.2126 * r + 0.7152 * g + 0.0722 * b) * 100

def check_palette_suitability(palette):
    """
    Checks if a palette is suitable for monochromatic vision.
    Returns a boolean and the list of calculated luminances.
    A palette is suitable if luminances are unique and separated by at least 5 points.
    """
    if len(palette) < 2:
        return True, []
    
    luminances = [calculate_luminance(hex_to_rgb_normalized(c)) for c in palette]
    
    # Check if rounded luminance values are unique
    rounded_lums = [round(lum) for lum in luminances]
    if len(set(rounded_lums)) != len(luminances):
        return False, luminances

    # Check if the minimum step between sorted luminances is large enough
    sorted_lums = sorted(luminances)
    min_step = min(sorted_lums[i+1] - sorted_lums[i] for i in range(len(sorted_lums) - 1))
    
    return min_step >= 5.0, luminances

# --- Color Palettes from the R code ---
# Plot 1 & 6 use the ggplot2 default palette, which has constant luminance.
pal_ggplot = ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']
# Plot 2 uses a diverging palette, which has varying luminance.
pal1 = ['#08519c', '#4292c6', '#f7f7f7', '#ef8a62', '#b2182b']
# Plot 3 uses HSLuv colors with constant Lightness (L=60) but varying Saturation.
pal2 = ['#999999', '#a29d88', '#82a884', '#3aaead', '#5e9bf3']
# Plot 4 uses HSLuv colors with constant Lightness (L=60) and Saturation (S=10).
pal3 = ['#a09797', '#a09a92', '#999d96', '#959d9c', '#9999a3']
# Plot 5 uses HSLuv with Lightness values sampled from a set with duplicates ({20,20,50,70,90}).
pal4 = ['#3c3838', '#868a7e', '#abbba7', '#e1e1e9', '#3b383d']

plots = {
    1: pal_ggplot,
    2: pal1,
    3: pal2,
    4: pal3,
    5: pal4,
    6: pal_ggplot
}

suitable_plots = []

print("Analyzing plot palettes for monochromatic visibility...\n")

for plot_num, palette in plots.items():
    is_suitable, luminances = check_palette_suitability(palette)
    
    print(f"Plot {plot_num}:")
    # The calculated luminance for each color in the palette
    lum_strings = [f"{lum:.1f}" for lum in luminances]
    print(f"  Luminance values = {', '.join(lum_strings)}")
    
    if is_suitable:
        print("  Result: SUITABLE")
        suitable_plots.append(str(plot_num))
    else:
        print("  Result: Not suitable")
    print("-" * 20)

if not suitable_plots:
    final_answer = "none"
else:
    final_answer = ",".join(sorted(list(set(suitable_plots))))

print(f"\nFinal Answer: The plot number(s) using a suitable color palette are: {final_answer}")