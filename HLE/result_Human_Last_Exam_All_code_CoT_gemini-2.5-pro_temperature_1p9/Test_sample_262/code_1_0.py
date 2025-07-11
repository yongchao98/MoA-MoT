import colorspacious
import hsluv
import numpy as np

def get_lightness(hex_color):
    """Converts a hex color string to its CIELab L* (lightness) value."""
    try:
        rgb_normalized = colorspacious.cspace_convert(hex_color, "hex", "sRGB1")
        lab = colorspacious.cspace_convert(rgb_normalized, "sRGB1", "CIELab")
        return lab[0]
    except Exception:
        return -1 # Should not happen with valid hex

def analyze_palettes():
    """
    Analyzes the color palettes from the plots to determine their
    suitability for monochromatic vision.
    """
    # Palette for Plot 1 and 6 (ggplot2 default / scales::hue_pal())
    # This palette has colors with the same luminance/lightness.
    # L* value for all colors is ~70
    pal_1_6 = [hsluv.hcl_to_hex((h, 50, 70)) for h in np.linspace(15, 375, 6)[:-1]]

    # Palette for Plot 2 (pals::ocean.balance(5))
    # These are the actual hex codes produced by the R package
    pal_2 = ["#1E4258", "#539E97", "#FFFFFF", "#F38548", "#7F3222"]

    # Palette for Plot 3
    # L is constant at 60, so L* will also be constant
    pal_3_h = np.arange(0, 300, 60)[:5]
    pal_3_s = np.array([h / 3 for h in pal_3_h])
    pal_3 = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip(pal_3_h, pal_3_s)]

    # Palette for Plot 4
    # L is constant at 60, so L* will also be constant
    pal_4_h = np.arange(0, 300, 60)[:5]
    pal_4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in pal_4_h]

    # Palette for Plot 5
    # The lightness values `l` are a permutation of (20, 50, 70, 20, 90)
    # This means two colors will have the same lightness value of 20.
    pal_5_h = np.arange(0, 289, 72)
    pal_5_l = [20.0, 50.0, 70.0, 20.0, 90.0]
    pal_5 = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(pal_5_h, pal_5_l)]
    
    palettes = {
        1: pal_1_6,
        2: pal_2,
        3: pal_3,
        4: pal_4,
        5: pal_5,
        6: pal_1_6
    }
    
    interpretable_plots = []
    print("Lightness (L*) analysis for each plot's palette:")
    print("A suitable palette must have distinct L* values for all colors.\n")

    for plot_num, palette in palettes.items():
        lightness_values = [round(get_lightness(c), 1) for c in palette]
        # Check if all L* values are unique
        is_interpretable = len(lightness_values) == len(set(lightness_values))
        
        print(f"Plot {plot_num}:")
        print(f"  Lightness values: {lightness_values}")
        
        if is_interpretable:
            print("  Result: Interpretable. All lightness values are unique.")
            interpretable_plots.append(str(plot_num))
        else:
            print("  Result: NOT interpretable. Contains non-unique lightness values.")
        print("-" * 25)

    print("\n--Conclusion--")
    if not interpretable_plots:
        print("None of the plots use a color palette interpretable for someone with full monochromatic vision.")
        final_answer = "none"
    else:
        final_answer = ",".join(sorted(interpretable_plots))
        print(f"The only plot using a suitable color palette is: {final_answer}")
    
if __name__ == '__main__':
    # The colorspacious and hsluv libraries are required.
    # You can install them with: pip install colorspacious hsluv
    try:
        analyze_palettes()
    except ImportError:
        print("Please install required libraries: pip install colorspacious hsluv")
    except Exception as e:
        print(f"An error occurred: {e}")