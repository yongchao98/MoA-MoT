def hex_to_luminance(hex_color):
    """Converts a hex color string to its perceptual luminance."""
    hex_color = hex_color.lstrip('#')
    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    # Standard formula for luminance
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return luminance

def analyze_palettes():
    """
    Analyzes color palettes to determine their suitability for monochromatic vision.
    A palette is suitable if its colors have distinct and well-separated luminance values.
    """
    # From R's: pals::ocean.balance(5)
    # This is a diverging palette, which is designed to vary in lightness.
    pal2_colors = ["#436384", "#98AABF", "#F2EBE5", "#C29270", "#8B4513"]

    # From R's: scales::hue_pal()(5)
    # The default ggplot palette has similar luminance values by design.
    pal6_colors = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

    print("Analyzing Plot 2 Palette (pals::ocean.balance):")
    luminances_2 = [hex_to_luminance(c) for c in pal2_colors]
    print(f"  Colors: {pal2_colors}")
    print(f"  Luminance values: {[round(l) for l in luminances_2]}")
    print("  Sorted Luminance: ", sorted([round(l) for l in luminances_2]))
    print("  Analysis: The luminance values are distinct and well-spaced, ranging from dark (84) to very light (236). This palette would be interpretable.\n")

    print("Analyzing Plots 1 & 6 Palette (scales::hue_pal):")
    luminances_6 = [hex_to_luminance(c) for c in pal6_colors]
    print(f"  Colors: {pal6_colors}")
    print(f"  Luminance values: {[round(l) for l in luminances_6]}")
    print("  Sorted Luminance: ", sorted([round(l) for l in luminances_6]))
    print("  Analysis: The luminance values are clustered in a narrow mid-range (126-160). This would be very difficult to interpret for someone with monochromatic vision.\n")
    
    print("Analysis of other plots:")
    print("  Plot 3 & 4: Use constant lightness (L=60) in HSLuv space, so colors are indistinguishable in grayscale.")
    print("  Plot 5: Uses lightness values from {20, 50, 70, 20, 90}, which contains a duplicate. Two colors will have identical lightness (20), making them indistinguishable.")
    
    print("\nConclusion: Only Plot 2 uses a suitable color palette.")
    print("\nFinal Answer:")


analyze_palettes()
print("2")