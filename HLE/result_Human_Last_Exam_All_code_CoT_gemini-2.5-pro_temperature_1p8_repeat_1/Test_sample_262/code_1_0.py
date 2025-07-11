import colorspacious
import hsluv
import sys

# This script analyzes color palettes to determine their suitability for monochromatic vision.
# A person with monochromatic vision perceives the world in shades of gray, so for a
# palette to be effective, each color must have a unique lightness. We test this by
# converting each color to the CIELAB color space and checking its L* (lightness) value.

# Note: The R code in the prompt generates some palettes procedurally.
# This script uses Python equivalents to recreate those palettes.

def analyze_palette(plot_num, palette_hex_codes, description):
    """
    Analyzes a palette, prints the results, and returns True if suitable, False otherwise.
    
    A palette is suitable if all its colors have distinct lightness values.
    """
    print(f"--- Analyzing Plot {plot_num}: {description} ---")
    
    # Helper function to convert a hex color string to a CIELAB L* value.
    def get_lightness(hex_color):
        try:
            # Convert hex to an RGB triplet normalized to 0-1
            rgb_normalized = tuple(int(hex_color.lstrip('#')[i:i+2], 16) / 255.0 for i in (0, 2, 4))
            # Convert sRGB to CIELab and return the L* component
            lab = colorspacious.cspace_convert(rgb_normalized, "sRGB1", "CIELab")
            return lab[0]
        except Exception as e:
            print(f"Error processing color {hex_color}: {e}", file=sys.stderr)
            return None

    lightness_values = [get_lightness(c) for c in palette_hex_codes]
    
    # Filter out any Nones from conversion errors
    lightness_values = [l for l in lightness_values if l is not None]
    
    # Round lightness to one decimal place for robust comparison
    rounded_lightness = [round(l, 1) for l in lightness_values]
    
    print(f"Colors (HEX): {palette_hex_codes}")
    print(f"Corresponding CIELAB L* values: {rounded_lightness}")
    
    # Check for suitability: The number of unique lightness values must equal the number of colors.
    is_suitable = len(set(rounded_lightness)) == len(palette_hex_codes)
    
    if is_suitable:
        print("Result: SUITABLE. All lightness values are distinct.")
    else:
        print("Result: NOT SUITABLE. Two or more colors have the same or very similar lightness.")
    print("-" * 50)
    
    return is_suitable, plot_num

def main():
    print("Step-by-step analysis of each plot's color palette:\n")

    # --- Palette Definitions ---

    # Plots 1 & 6: `ggplot` default / `scales::hue_pal()`. Uses HCL color space with constant
    # luminance (l=65) and chroma, which results in colors with very similar lightness.
    # Generating 5 colors as used in the plot.
    pal_1_and_6 = [hsluv.hsluv_to_hex((h, 100, 65)) for h in [15.0, 88.0, 161.0, 234.0, 307.0]]

    # Plot 2: `pals::ocean.balance(5)`. A diverging palette known to be colorblind-friendly.
    pal_2 = ["#4B3A82", "#4988A9", "#E2E2E2", "#E29847", "#B43D2A"]

    # Plot 3: Custom HSLuv. `l=60` is constant, making it unsuitable by definition.
    hues_3 = [0, 60, 120, 180, 240] # Using first 5 values for 5 groups
    saturations_3 = [h/3.0 for h in hues_3]
    pal_3 = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip(hues_3, saturations_3)]

    # Plot 4: Custom HSLuv. `l=60` is constant, also unsuitable.
    pal_4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in [0, 60, 120, 180, 240]]

    # Plot 5: Custom HSLuv. Lightness is sampled from a set with a duplicate value (20).
    # This guarantees at least two colors will have the same lightness.
    hues_5 = [0, 72, 144, 216, 288]
    lightnesses_5 = [20, 50, 70, 20, 90]
    pal_5 = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(hues_5, lightnesses_5)]
    
    palettes = {
        1: (pal_1_and_6, "Default ggplot / scales::hue_pal()"),
        2: (pal_2, "pals::ocean.balance diverging palette"),
        3: (pal_3, "Custom HSLuv with constant lightness (L=60)"),
        4: (pal_4, "Custom HSLuv with constant lightness (L=60)"),
        5: (pal_5, "Custom HSLuv with duplicate lightness values"),
        6: (pal_1_and_6, "Default ggplot / scales::hue_pal()")
    }

    # --- Run Analysis ---
    results = []
    for num, (pal, desc) in palettes.items():
        results.append(analyze_palette(num, pal, desc))
    
    suitable_plots = [str(num) for is_suitable, num in results if is_suitable]

    final_answer = ",".join(suitable_plots) if suitable_plots else "none"

    print(f"\nFinal Answer: The plot number(s) interpretable under monochromatic vision are: {final_answer}")


if __name__ == '__main__':
    main()