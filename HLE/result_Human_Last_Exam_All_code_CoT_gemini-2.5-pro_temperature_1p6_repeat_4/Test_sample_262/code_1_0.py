import sys
import subprocess

# This script uses 'colorspacious' and 'hsluv' to analyze colors.
# This block will install them if they are not already in your environment.
try:
    import colorspacious
    import hsluv
except ImportError:
    print("Installing required Python libraries: colorspacious, hsluv...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "colorspacious", "hsluv"],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        import colorspacious
        import hsluv
        print("Libraries installed successfully.")
    except Exception as e:
        print(f"Error installing libraries: {e}")
        print("Please install them manually using: pip install colorspacious hsluv")
        sys.exit(1)

def hex_to_rgb_normalized(hex_color):
    """Converts a hex color string to an (R, G, B) tuple normalized to 0-1."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))

def get_cielab_lightness(hex_color):
    """Converts a hex color to its CIELab L* (perceptual lightness) value."""
    # The standard conversion path is sRGB (for screens) -> XYZ -> CIELab
    rgb_srgb1 = hex_to_rgb_normalized(hex_color)
    xyz = colorspacious.cspace_convert(rgb_srgb1, "sRGB1", "XYZ100")
    lab = colorspacious.cspace_convert(xyz, "XYZ100", "CIELab")
    return lab[0]

def analyze_palettes():
    """
    Analyzes the color palettes from the problem description to determine
    their suitability for monochromatic vision.
    """
    # --- Step 1: Define the Palettes based on the R code ---

    # Plot 1 & 6: Default ggplot2 palette from scales::hue_pal()(5)
    pal1_and_6 = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

    # Plot 2: Diverging palette from pals::ocean.balance(5)
    pal2_ocean = ["#1B7837", "#76B462", "#E8E8E8", "#98569D", "#4D004B"]

    # Plot 3: HSLuv with constant lightness (60), first 5 colors used
    hues_p3 = [i for i in range(0, 300, 60)]
    sats_p3 = [i / 3 for i in hues_p3]
    pal3_hsluv = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip(hues_p3, sats_p3)][:5]

    # Plot 4: HSLuv with constant lightness (60), first 5 colors used
    hues_p4 = [i for i in range(0, 300, 60)][:5]
    pal4_hsluv = [hsluv.hsluv_to_hex((h, 10, 60)) for h in hues_p4]

    # Plot 5: HSLuv with sampled lightness, which includes a duplicate value (20)
    hues_p5 = [i for i in range(0, 289, 72)]
    lightness_p5 = [20, 50, 70, 20, 90]
    pal5_hsluv = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(hues_p5, lightness_p5)]

    palettes = {
        1: pal1_and_6,
        2: pal2_ocean,
        3: pal3_hsluv,
        4: pal4_hsluv,
        5: pal5_hsluv,
        6: pal1_and_6,
    }

    # --- Step 2: Analyze each palette ---
    print("A plot is interpretable in monochrome if the colors in its palette")
    print("have distinct perceptual lightness (L*) values.\n")

    suitable_plots = []
    for plot_num in sorted(palettes.keys()):
        palette_colors = palettes[plot_num]
        lightness_values = [get_cielab_lightness(c) for c in palette_colors]
        # Rounding values helps compare them robustly.
        rounded_lightness = [round(l, 1) for l in lightness_values]
        
        # Check if the number of unique lightness values equals the number of colors.
        is_suitable = len(set(rounded_lightness)) == len(rounded_lightness)

        print(f"--- Analysis for Plot {plot_num} ---")
        print(f"Palette Colors: {palette_colors}")
        print(f"Perceptual Lightness (L*) values: {[f'{l:.1f}' for l in lightness_values]}")
        
        if is_suitable:
            print("Result: SUITABLE. All lightness values are unique.")
            suitable_plots.append(str(plot_num))
        else:
            print("Result: NOT SUITABLE. At least two colors have the same lightness.")
        print()
    
    # --- Step 3: Conclude and print the final answer ---
    if not suitable_plots:
        final_answer = "none"
    else:
        # Sort for consistent output, though there's only one correct answer here.
        final_answer = ",".join(sorted(suitable_plots))
        
    print(f"Final Conclusion: The only plot that uses a color palette fully interpretable")
    print(f"for someone with monochromatic vision is: {final_answer}")

if __name__ == "__main__":
    analyze_palettes()
<<<2>>>