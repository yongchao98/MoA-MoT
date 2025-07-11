import sys
try:
    from colorspacious import cspace_convert
    import hsluv
except ImportError:
    print("Please install the required libraries by running: pip install colorspacious hsluv")
    sys.exit(1)

def hex_to_srgb1(hex_color):
    """Converts a hex color string to an sRGB1 tuple (values 0-1)."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))

def get_grayscale_value(srgb1_color):
    """
    Simulates monochromatic vision to get a single grayscale lightness value for a color.
    We convert to the CAM02-UCS color space and take the J' (Lightness) component.
    """
    return cspace_convert(srgb1_color, "sRGB1", "CAM02-UCS")['J\'']

def analyze_palette(name, hex_colors, explanation=""):
    """Analyzes a palette for monochromatic visibility and prints the results."""
    print(f"--- Analyzing {name} ---")
    if explanation:
        print(explanation)
    
    lightness_values = [round(get_grayscale_value(hex_to_srgb1(h)), 2) for h in hex_colors]
    
    # Check for uniqueness of lightness values
    is_interpretable = len(set(lightness_values)) == len(hex_colors)
    
    print(f"Colors: {hex_colors}")
    print(f"Perceived lightness values: {lightness_values}")
    
    if is_interpretable:
        print("Result: This palette IS SUITABLE for monochromatic vision.")
        return True
    else:
        print("Result: This palette IS NOT SUITABLE as some colors have the same lightness.")
        return False

def main():
    # Palette definitions from the R code
    palettes = {
        "Plot 1 & 6 (ggplot default)": {
            "colors": ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"],
            "explanation": "These palettes use colors of varying hue but similar lightness, which is typical for qualitative palettes."
        },
        "Plot 2 (pals::ocean.balance)": {
            "colors": ["#0A42A4", "#7199D5", "#FFFFFF", "#F38C6C", "#D2222E"],
            "explanation": "This is a divergent palette. Such palettes are often designed to vary in lightness from dark to light to dark, making them readable in grayscale."
        },
        "Plot 3 (pal2 from hsluv)": {
            "colors": [hsluv.hsluv_to_hex([h, s, 60]) for h, s in zip([0, 60, 120, 180, 240], [0, 20, 40, 60, 80])],
            "explanation": "This palette was generated with a constant HSLuv Lightness (l=60). HSLuv is designed so that colors with the same 'l' value have the same perceived lightness."
        },
        "Plot 4 (pal3 from hsluv)": {
            "colors": [hsluv.hsluv_to_hex([h, 10, 60]) for h in [0, 60, 120, 180, 240]],
            "explanation": "Like Plot 3, this palette was also generated with a constant HSLuv Lightness (l=60), making the colors indistinguishable in grayscale."
        },
        "Plot 5 (pal4 from hsluv)": {
            "colors": [hsluv.hsluv_to_hex([h, 10, l]) for h, l in zip([0, 72, 144, 216, 288], [20, 50, 70, 20, 90])],
            "explanation": "This palette was generated using the HSLuv lightness values [20, 50, 70, 20, 90]. Since the lightness value '20' is used twice, two of the colors will have the same perceived lightness."
        }
    }
    
    suitable_plots = []
    
    for name, data in palettes.items():
        # The key for Plot 1 and 6 is the same, so we add both if it's suitable.
        plot_numbers = [int(s) for s in name.split() if s.isdigit()]

        print("") # Newline for readability
        if analyze_palette(name, data['colors'], data['explanation']):
            suitable_plots.extend(plot_numbers)
    
    print("\n--- Conclusion ---")
    if suitable_plots:
        result = ",".join(map(str, sorted(suitable_plots)))
        print(f"The plots that use a color palette interpretable for someone with full monochromatic vision are: {result}")
    else:
        result = "none"
        print("None of the plots use a color palette interpretable for someone with full monochromatic vision.")
    
    # Final answer in requested format
    print(f"\n<<<{result}>>>")


if __name__ == "__main__":
    main()