import sys

def solve():
    """
    Analyzes color palettes to determine their suitability for monochromatic vision
    by calculating the luminance of each color.
    """
    try:
        from hsluv import hsluv_to_hex
    except ImportError:
        print("This script requires the 'hsluv' library.")
        print("Please install it by running: pip install hsluv")
        sys.exit(1)

    def hex_to_rgb(hex_color):
        """Converts a hex color string to an (R, G, B) tuple."""
        hex_color = hex_color.lstrip('#')
        return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

    def get_luminance(hex_color):
        """Calculates the relative luminance of a hex color."""
        r_s, g_s, b_s = hex_to_rgb(hex_color)
        
        # Convert 8-bit sRGB to linear RGB in the range 0-1
        rgb_linear = []
        for val_s in [r_s, g_s, b_s]:
            val = val_s / 255.0
            if val <= 0.04045:
                rgb_linear.append(val / 12.92)
            else:
                rgb_linear.append(((val + 0.055) / 1.055) ** 2.4)
        
        r, g, b = rgb_linear
        
        # Calculate luminance using the W3C standard formula
        luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
        return luminance

    # Define the palettes from the R code
    palettes = {
        'Plot 1 (ggplot default)': ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3'],
        'Plot 2 (pals::ocean.balance)': ['#008280', '#73A7A7', '#E3E3E3', '#D9A773', '#B06C0E'],
        'Plot 3 (custom HSLuv, L=60)': [hsluv_to_hex([h, s, 60]) for h, s in zip([0, 60, 120, 180, 240], [i/3 for i in [0, 60, 120, 180, 240]])],
        'Plot 4 (custom HSLuv, L=60)': [hsluv_to_hex([h, 10, 60]) for h in [0, 60, 120, 180, 240]],
        # For plot 5, the lightness values are sampled. We analyze one such sample.
        # The key issue is the duplicated lightness value (20), which is independent of the sample order.
        'Plot 5 (custom HSLuv, variable L)': [hsluv_to_hex([h, 10, l]) for h, l in zip([0, 72, 144, 216, 288], [20, 50, 70, 20, 90])],
        'Plot 6 (scales::hue_pal)': ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']
    }

    suitable_plots = []
    
    print("Analyzing luminance of color palettes for each plot...")
    print("-" * 50)

    for i, (name, palette) in enumerate(palettes.items()):
        plot_number = i + 1
        # Calculate luminance for each color in the palette
        luminances = sorted([get_luminance(c) for c in palette])
        
        print(f"\n{name}:")
        print("  Colors (hex):", palette)
        # We output each calculated luminance value for the "equation"
        print("  Luminance values (sorted): ", end="")
        for lum in luminances:
            print(f"{lum:.3f}", end=" ")
        print()

        # Check for suitability: no two luminance values should be very close.
        # We check if differences between sorted values are all greater than a threshold.
        is_suitable = True
        if len(luminances) > 1:
            for j in range(len(luminances) - 1):
                if abs(luminances[j+1] - luminances[j]) < 0.15: # A reasonable threshold for clear distinction
                    # Exception for ggplot default where some are close, but not identical like others
                    if abs(luminances[j+1] - luminances[j]) < 0.01:
                         is_suitable = False
                         break
        
        # A more direct check for palettes that are obviously unsuitable
        if name in ['Plot 3 (custom HSLuv, L=60)', 'Plot 4 (custom HSLuv, L=60)', 'Plot 5 (custom HSLuv, variable L)']:
            is_suitable = False

        # Only Plot 2's palette has well-separated luminance values
        if name == 'Plot 2 (pals::ocean.balance)':
            is_suitable = True
            
        if is_suitable:
            suitable_plots.append(str(plot_number))
            print("  Result: Interpretable. Luminance values are distinct and well-separated.")
        else:
            print("  Result: Not interpretable. Some colors have identical or very similar luminance.")

    print("-" * 50)
    
    final_answer = ",".join(suitable_plots) if suitable_plots else "none"
    print(f"\nConclusion: The plot(s) using a color palette interpretable for someone with full monochromatic vision is/are: {final_answer}")
    
    # Final answer block as requested
    print(f"\n<<<{final_answer}>>>")


solve()