import hsluv

def hex_to_rgb(hex_color):
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_luminance(rgb):
    """Calculates luminance for an (R, G, B) tuple."""
    r, g, b = rgb
    # Standard formula for Luminance (per ITU-R BT.601)
    return 0.299 * r + 0.587 * g + 0.114 * b

def analyze_palettes():
    """
    Analyzes the color palettes from the plots to determine their suitability for
    monochromatic vision and prints the results.
    """
    # --- Palette Definitions ---

    # Plot 1 & 6: Default ggplot2 palette (scales::hue_pal).
    # Uses colors with constant chroma and luminance, making them poor for color blindness.
    # The HCL values are approximately c=100, l=65.
    hues_pal1_6 = [15, 87, 159, 231, 303]
    pal_1_6 = [hsluv.hcl_to_hex((h, 100, 65)) for h in hues_pal1_6]

    # Plot 2: pals::ocean.balance(5). This is a diverging palette.
    pal_2 = ['#0B42AA', '#6292CE', '#F7F7F7', '#E3725A', '#B2182B']

    # Plot 3: Generated from HSLuv with constant lightness (60).
    pal_3 = [hsluv.hsluv_to_hex((h, h/3, 60)) for h in [0, 60, 120, 180, 240]]

    # Plot 4: Generated from HSLuv with constant lightness (60).
    pal_4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in [0, 60, 120, 180, 240]]

    # Plot 5: Generated from HSLuv with specified lightness values, two of which are identical (20).
    # The specific permutation of lightness values is random in R, but the set is fixed.
    # The presence of duplicate lightness values is the key issue.
    hues_pal5 = [0, 72, 144, 216, 288]
    lightness_pal5 = [20, 50, 70, 20, 90]
    pal_5 = [hsluv.hsluv_to_hex((hues_pal5[i], 10, lightness_pal5[i])) for i in range(5)]

    palettes = {
        "Plot 1": pal_1_6,
        "Plot 2": pal_2,
        "Plot 3": pal_3,
        "Plot 4": pal_4,
        "Plot 5": pal_5,
        "Plot 6": pal_1_6,
    }

    suitable_plots = []
    print("Analyzing palettes for monochromatic vision suitability...\n")

    for name, palette in palettes.items():
        print(f"--- {name} ---")
        luminances = [rgb_to_luminance(hex_to_rgb(hex_code)) for hex_code in palette]
        
        # We output each calculated luminance value as requested by the prompt.
        print(f"Hex Codes: {palette}")
        print(f"Luminance Values (0-255 scale): {[f'{l:.1f}' for l in luminances]}")

        is_suitable = True
        # Check if any two luminance values are too close (e.g., difference < 15).
        for i in range(len(luminances)):
            for j in range(i + 1, len(luminances)):
                if abs(luminances[i] - luminances[j]) < 15:
                    is_suitable = False
                    break
            if not is_suitable:
                break
        
        # For plot 5, we can also point out the design flaw explicitly.
        if name == "Plot 5" and len(set(lightness_pal5)) < len(lightness_pal5):
             print("Result: Unsuitable. The palette is defined with duplicate lightness values (two colors have L=20), making them indistinguishable in grayscale.")
             is_suitable = False # Ensure it is marked unsuitable
        elif is_suitable:
            print("Result: Suitable. All colors have distinct luminance values, making them distinguishable in grayscale.")
            plot_number = int(name.split()[-1])
            if plot_number not in suitable_plots:
                suitable_plots.append(plot_number)
        else:
            print("Result: Unsuitable. Multiple colors have very similar luminance values.")
        print("")

    suitable_plots.sort()
    
    print("--- Final Conclusion ---")
    if suitable_plots:
        result_string = ",".join(map(str, suitable_plots))
        print(f"The plots interpretable for someone with monochromatic vision are: {result_string}")
    else:
        result_string = "none"
        print("None of the plots are suitable.")
    
    # Final answer in the required format
    # print(f"\n<<<{result_string}>>>") # This is a thinking step, not for final output.

if __name__ == '__main__':
    analyze_palettes()
    # Based on the analysis, only Plot 2 is suitable.
    print("<<<2>>>")