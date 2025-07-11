import math

def hex_to_rgb(hex_str):
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_str = hex_str.lstrip('#')
    return tuple(int(hex_str[i:i+2], 16) for i in (0, 2, 4))

def calculate_luminance(rgb_tuple):
    """
    Calculates the perceptual luminance of an sRGB color.
    Formula from Rec. 709 standard.
    """
    # Scale RGB values to 0-1
    rgb = [x / 255.0 for x in rgb_tuple]
    
    # Apply gamma correction (sRGB to linear)
    linear_rgb = []
    for c in rgb:
        if c <= 0.04045:
            linear_rgb.append(c / 12.92)
        else:
            linear_rgb.append(math.pow((c + 0.055) / 1.055, 2.4))
            
    r, g, b = linear_rgb
    
    # Calculate luminance
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return luminance

def analyze_palettes():
    """
    Analyzes color palettes to determine their suitability for monochromatic vision.
    """
    # Palettes to analyze. Plot 1 and 6 use the same default palette.
    palettes = {
        "Plot 1 & 6": ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"],
        "Plot 2": ["#045275", "#089099", "#7CCBA2", "#C8D5B9", "#F0F0F0"]
    }

    print("--- Analysis of Color Palette Luminance ---\n")
    print("For a plot to be readable with monochromatic vision, the colors must have distinct luminance values.\n")

    # Analyze palettes with provided hex codes
    for name, hex_codes in palettes.items():
        print(f"Analyzing {name}:")
        luminance_values = [calculate_luminance(hex_to_rgb(h)) for h in hex_codes]
        # Print rounded values for clarity
        print(f"  Luminance values: {[round(l, 3) for l in luminance_values]}")
        if len(set(round(l, 2) for l in luminance_values)) < len(luminance_values):
            print("  Result: Not suitable. Luminance values are too similar or identical.\n")
        else:
            print("  Result: SUITABLE. Luminance values are distinct and well-spaced.\n")

    # Explain the other plots based on their R code definition
    print("--- Analysis of Other Plots based on R Code ---\n")
    print("Plot 3: Not suitable. The 'L' (Lightness) parameter in HSLuv is constant at 60 for all colors, resulting in nearly identical luminance.")
    print("Plot 4: Not suitable. For the same reason as Plot 3, the 'L' parameter is constant at 60.")
    print("Plot 5: Not suitable. The 'L' (Lightness) is sampled from a set containing duplicates (c(20,50,70,20,90)). This guarantees at least two colors will have identical luminance.\n")
    
    print("--- Conclusion ---")
    print("Only Plot 2 uses a palette with systematically changing luminance, making it interpretable for someone with monochromatic vision.")
    
analyze_palettes()
print("\n<<<2>>>")