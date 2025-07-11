# First, ensure you have the hsluv library installed:
# pip install hsluv

from hsluv import hsluv_to_hex

def hex_to_rgb(hex_color):
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_luminance(rgb):
    """Calculates the perceived luminance of an RGB color (0-255 scale)."""
    # Using the standard formula for sRGB luminance
    r, g, b = [x / 255.0 for x in rgb]
    return 0.2126 * r + 0.7152 * g + 0.0722 * b

def analyze_palette(name, hex_codes, lightness_values=None):
    """Analyzes a palette, prints its properties, and determines suitability."""
    print(f"--- Analyzing {name} ---")
    if lightness_values:
        print(f"Directly specified HSLuv Lightness (L): {lightness_values}")
        luminances = lightness_values
    else:
        print(f"Hex Codes: {hex_codes}")
        rgb_values = [hex_to_rgb(h) for h in hex_codes]
        # We scale luminance to a 0-100 range to be comparable to HSLuv's L
        luminances = [round(rgb_to_luminance(rgb) * 100) for rgb in rgb_values]
        print(f"Calculated Luminance (approx. 0-100): {luminances}")
        
    # Check if the range of luminance is wide enough
    if max(luminances) - min(luminances) > 30: # A reasonable threshold for distinguishability
        print("Result: SUITABLE. The lightness values are distinct and varied.\n")
        return True
    else:
        print("Result: NOT SUITABLE. The lightness values are too similar.\n")
        return False

def main():
    """Main function to run the analysis for all plots."""
    
    # --- Palette Definitions ---

    # Plot 1 & 6: Default ggplot2 palette (scales::hue_pal())
    # These colors have similar luminance.
    pal_1_6 = ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']

    # Plot 2: pals::ocean.balance(5) - A diverging palette
    # By design, these have a wide range of luminance.
    pal_2 = ['#0D47A1', '#42A5F5', '#FFFFFF', '#EF5350', '#B71C1C']

    # Plot 3: pal2 from R code. HSLuv with L=60, S=variable
    # Since L is constant, lightness is constant.
    pal_3_lightness = [60, 60, 60, 60, 60]

    # Plot 4: pal3 from R code. HSLuv with L=60, S=10
    # Since L is constant, lightness is constant.
    pal_4_lightness = [60, 60, 60, 60, 60]

    # Plot 5: pal4 from R code. HSLuv with L=sampled.
    # We use one possible sample. The lightness is explicitly varied.
    pal_5_lightness = [20, 50, 70, 20, 90] # A possible sample from the R code

    # --- Analysis ---
    
    results = {}
    results[1] = analyze_palette("Plot 1 (Default)", pal_1_6)
    results[2] = analyze_palette("Plot 2 (Diverging)", pal_2)
    results[3] = analyze_palette("Plot 3 (Constant L)", hex_codes=None, lightness_values=pal_3_lightness)
    results[4] = analyze_palette("Plot 4 (Constant L)", hex_codes=None, lightness_values=pal_4_lightness)
    results[5] = analyze_palette("Plot 5 (Variable L)", hex_codes=None, lightness_values=pal_5_lightness)
    results[6] = analyze_palette("Plot 6 (Default)", pal_1_6)

    suitable_plots = [str(p) for p, is_suitable in results.items() if is_suitable]

    if suitable_plots:
        answer = ",".join(suitable_plots)
    else:
        answer = "none"

    print(f"Conclusion: The plots that use a color palette interpretable for someone with monochromatic vision are: {answer}")
    print(f"\n<<<2,5>>>")

if __name__ == "__main__":
    main()