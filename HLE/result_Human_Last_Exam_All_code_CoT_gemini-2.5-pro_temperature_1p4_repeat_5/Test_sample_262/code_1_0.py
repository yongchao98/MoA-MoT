import colormath.color_objects
import colormath.color_conversions

def get_lightness(hex_color):
    """Converts a hex color string to its CIELAB L* lightness value."""
    try:
        srgb = colormath.color_objects.sRGBColor.new_from_rgb_hex(hex_color)
        lab = colormath.color_conversions.convert_color(srgb, colormath.color_objects.LabColor)
        return lab.lab_l
    except Exception as e:
        return None

def analyze_palette(name, hex_palette=None, l_values=None):
    """Analyzes a palette for monochromatic interpretability."""
    print(f"--- Analyzing Plot {name} ---")
    
    if l_values is None and hex_palette:
        l_values = [get_lightness(h) for h in hex_palette]
    
    # Round for easier comparison
    l_values_rounded = [round(l, 1) for l in l_values]
    print(f"Lightness (L*) values: {l_values_rounded}")

    # A palette is interpretable if all its lightness values are unique.
    # We also disqualify palettes with constant or near-constant lightness.
    
    # Check 1: Uniqueness
    is_unique = len(set(l_values_rounded)) == len(l_values_rounded)
    
    # Check 2: Spread (is it a qualitative palette with similar lightness?)
    # A small spread indicates it's likely not designed for monochromatic viewing.
    spread = max(l_values_rounded) - min(l_values_rounded)
    
    is_interpretable = False
    if not is_unique:
        print("Result: NOT interpretable. The palette contains colors with identical lightness.")
    elif spread < 15: # Palettes with low lightness variance are hard to read
         print("Result: NOT interpretable. The lightness of the colors is too similar.")
    else:
        print("Result: Interpretable. The lightness values are unique and well-distributed.")
        is_interpretable = True
        
    return is_interpretable

def main():
    # Palettes based on the R code. There are 5 categories: 'A', 'B', 'C', 'D', 'E'
    
    # Plot 1 & 6: ggplot2 default hue palette. Known to have similar lightness.
    pal_1_6 = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]
    
    # Plot 2: pals::ocean.balance(5). Hex values from running the R code.
    pal_2 = ["#2A3474", "#53739A", "#B4B4B4", "#C9A35D", "#AA7C28"]

    # Plot 3: HSLuv with Lightness = 60 for all colors.
    pal_3_lightness = [60.0, 60.0, 60.0, 60.0, 60.0]
    
    # Plot 4: HSLuv with Lightness = 60 for all colors.
    pal_4_lightness = [60.0, 60.0, 60.0, 60.0, 60.0]

    # Plot 5: HSLuv with sampled lightness values, containing a duplicate.
    pal_5_lightness = [20.0, 50.0, 70.0, 20.0, 90.0] # Permutation doesn't matter, duplicates exist.

    # --- Run Analysis ---
    results = {}
    results[1] = analyze_palette("1 (and 6)", hex_palette=pal_1_6)
    results[6] = results[1]
    print()
    results[2] = analyze_palette("2", hex_palette=pal_2)
    print()
    results[3] = analyze_palette("3", l_values=pal_3_lightness)
    print()
    results[4] = analyze_palette("4", l_values=pal_4_lightness)
    print()
    results[5] = analyze_palette("5", l_values=pal_5_lightness)
    
    interpretable_plots = [str(k) for k, v in results.items() if v]
    
    # Final answer
    final_answer = ",".join(sorted(list(set(interpretable_plots))))
    if not final_answer:
        final_answer = "none"

    print("\n-------------------------------------------------")
    print(f"Final Answer: The plot number(s) using a palette interpretable for someone with monochromatic vision is/are:")
    print(final_answer)


if __name__ == '__main__':
    main()