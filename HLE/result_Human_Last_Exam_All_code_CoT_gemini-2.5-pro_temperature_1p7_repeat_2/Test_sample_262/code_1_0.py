import matplotlib.colors as mcolors
import hsluv
import numpy as np

def srgb_to_linear(c_srgb):
    """Converts a single sRGB channel value (0-1) to a linear value."""
    if c_srgb <= 0.04045:
        return c_srgb / 12.92
    else:
        return ((c_srgb + 0.055) / 1.055) ** 2.4

def get_perceptual_luminance(hex_color):
    """
    Calculates the perceptual luminance for a given hex color code.
    This uses the standard formula after converting sRGB to linear RGB.
    """
    # Convert hex to an sRGB tuple (R, G, B) with values from 0-1
    rgb_srgb = mcolors.to_rgb(hex_color)
    
    # Convert sRGB to linear RGB
    rgb_linear = [srgb_to_linear(c) for c in rgb_srgb]
    
    # Calculate luminance using the standard weights
    luminance = (0.2126 * rgb_linear[0] + 
                 0.7152 * rgb_linear[1] + 
                 0.0722 * rgb_linear[2])
    return luminance

def analyze_palette(plot_num, palette_name, hex_colors):
    """Analyzes a palette and prints the results."""
    print(f"--- Analyzing Plot {plot_num}: {palette_name} ---")
    
    luminances = [get_perceptual_luminance(c) for c in hex_colors]
    
    # Round luminances for easier comparison
    rounded_lums = [round(l, 2) for l in luminances]
    
    # Check if there are duplicate luminance values
    is_interpretable = len(set(rounded_lums)) == len(hex_colors)
    
    print(f"Colors (Hex): {', '.join(hex_colors)}")
    print(f"Luminance values: {[f'{l:.3f}' for l in luminances]}")
    
    if not is_interpretable:
        print("Conclusion: NOT interpretable. Some colors have identical or very similar luminance.")
    else:
        # Also check if the values are well-spaced
        if np.std(luminances) < 0.1:
            print("Conclusion: Barely interpretable. Luminance values are unique but very close together, making them hard to distinguish.")
        else:
            print("Conclusion: Interpretable. Luminance values are distinct and well-separated.")
    print("-" * (18 + len(palette_name)))


def solve():
    """
    Solves the problem by defining, analyzing, and reporting on each palette.
    """
    # Palette definitions based on the R code
    
    # Plot 1 & 6: Default ggplot2 palette (scales::hue_pal()(5))
    # This generates hues evenly spaced around the color wheel at a constant lightness and chroma.
    pal1_6 = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

    # Plot 2: pals::ocean.balance(5) - A diverging palette
    pal2 = ["#0D4A7F", "#6895A9", "#F2F2F2", "#F2AC6C", "#BF512A"]

    # Plot 3: HSLuv with varying Hue and Saturation, constant Lightness (l=60)
    # The data 'g' has 5 levels, so ggplot uses the first 5 generated colors.
    hues_pal3 = [0, 60, 120, 180, 240]
    sats_pal3 = [0, 20, 40, 60, 80] # s = i/3 from the R code
    pal3 = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip(hues_pal3, sats_pal3)]

    # Plot 4: HSLuv with varying Hue, constant Saturation (s=10) and Lightness (l=60)
    # Colors with the same HSLuv Lightness 'L' are designed to be perceptually isoluminant (same brightness).
    hues_pal4 = [0, 60, 120, 180, 240]
    pal4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in hues_pal4]

    # Plot 5: HSLuv with varying Hue and sampled Lightness.
    # The lightness values are sampled from a list that contains a duplicate (20).
    hues_pal5 = [0, 72, 144, 216, 288]
    lightness_pal5 = [20, 50, 70, 20, 90] # One possible sample order
    pal5 = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(hues_pal5, lightness_pal5)]
    
    palettes = {
        1: ("Default ggplot2 palette", pal1_6),
        2: ("pals::ocean.balance", pal2),
        3: ("HSLuv (varying S, L=60)", pal3),
        4: ("HSLuv (S=10, L=60)", pal4),
        5: ("HSLuv (sampled L with duplicate)", pal5),
        6: ("Default ggplot2 palette", pal1_6)
    }
    
    interpretable_plots = []
    for plot_num, (name, colors) in palettes.items():
        analyze_palette(plot_num, name, colors)
        
        # Determine final answer based on strict criteria
        luminances = [get_perceptual_luminance(c) for c in colors]
        rounded_lums = [round(l, 2) for l in luminances]
        if len(set(rounded_lums)) == len(colors) and np.std(luminances) > 0.1:
            interpretable_plots.append(str(plot_num))
            
    print("\nBased on the analysis, only palettes with clearly distinct and well-separated luminance values are considered interpretable.")
    print("The plot number that meets these criteria is:")
    print(','.join(interpretable_plots) if interpretable_plots else "none")


solve()
