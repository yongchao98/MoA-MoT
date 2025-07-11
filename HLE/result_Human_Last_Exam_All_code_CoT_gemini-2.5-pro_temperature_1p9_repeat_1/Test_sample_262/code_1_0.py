import hsluv
import numpy as np

def hex_to_rgb(hex_str):
    """Converts a hex color string to an RGB tuple in [0, 1] range."""
    hex_str = hex_str.lstrip('#')
    return tuple(int(hex_str[i:i+2], 16) / 255.0 for i in (0, 2, 4))

def calculate_and_print_luminance(rgb_tuple, color_hex):
    """Calculates relative luminance and prints the formula and result."""
    r, g, b = rgb_tuple
    # Standard formula for relative luminance (Y)
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    print(f"Color {color_hex}: Luminance = 0.2126 * {r:.2f} + 0.7152 * {g:.2f} + 0.0722 * {b:.2f} = {luminance:.4f}")
    return luminance

def analyze_plots():
    """
    Analyzes the color palettes from the plots to determine which are suitable
    for viewers with monochromatic vision.
    """
    print("Analyzing plot palettes for monochromatic vision suitability.\n")
    suitable_plots = []

    # --- Analysis for Plot 1 & 6 (Default Hue-based Palette) ---
    print("--- Analysis: Plot 1 & 6 ---")
    print("Palette based on evenly spaced hues, which often have similar luminance.")
    pal_1_6 = ['#F8766D', '#7CAE00', '#00BFC4', '#C77CFF', '#F564E3']
    lums_1_6 = [calculate_and_print_luminance(hex_to_rgb(c), c) for c in pal_1_6]
    range_1_6 = max(lums_1_6) - min(lums_1_6)
    print(f"Result: NOT SUITABLE. The luminance range is very low ({range_1_6:.4f}).\n")

    # --- Analysis for Plot 2 (Diverging Palette) ---
    print("--- Analysis: Plot 2 ---")
    print("Diverging palette, designed to vary in lightness.")
    pal_2 = ['#0B4294', '#74A7C2', '#F7F7F7', '#E4A07B', '#B4532A']
    lums_2 = [calculate_and_print_luminance(hex_to_rgb(c), c) for c in pal_2]
    # Check if luminances are distinct (checking rounded values is a robust way)
    if len(set([round(l, 2) for l in lums_2])) == len(lums_2):
        print("Result: SUITABLE. Luminances are distinct and cover a wide range.\n")
        suitable_plots.append(2)
    else:
        print("Result: NOT SUITABLE.\n")
    
    # --- Analysis for Plot 3 (Constant Lightness, Varying Saturation) ---
    print("--- Analysis: Plot 3 ---")
    print("Palette with constant HSLuv lightness (L=60).")
    hues_3 = [0, 60, 120, 180, 240]
    sats_3 = [h / 3 for h in hues_3]
    pal_3 = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip(hues_3, sats_3)]
    lums_3 = [calculate_and_print_luminance(hex_to_rgb(c), c) for c in pal_3]
    range_3 = max(lums_3) - min(lums_3)
    print(f"Result: NOT SUITABLE. Constant HSLuv lightness results in near-identical luminance ({range_3:.4f} range).\n")
    
    # --- Analysis for Plot 4 (Constant Lightness and Saturation) ---
    print("--- Analysis: Plot 4 ---")
    print("Palette with constant HSLuv lightness (L=60) and saturation (S=10).")
    hues_4 = [0, 60, 120, 180, 240]
    pal_4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in hues_4]
    lums_4 = [calculate_and_print_luminance(hex_to_rgb(c), c) for c in pal_4]
    range_4 = max(lums_4) - min(lums_4)
    print(f"Result: NOT SUITABLE. Constant HSLuv lightness results in near-identical luminance ({range_4:.4f} range).\n")

    # --- Analysis for Plot 5 (Varying Lightness with Duplicates) ---
    print("--- Analysis: Plot 5 ---")
    print("Palette with HSLuv lightness values including duplicates.")
    lightness_values = [20, 50, 70, 20, 90]
    hues_5 = [0, 72, 144, 216, 288]
    print(f"Using HSLuv lightness values: {lightness_values}")
    pal_5 = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(hues_5, lightness_values)]
    lums_5 = [calculate_and_print_luminance(hex_to_rgb(c), c) for c in pal_5]
    if len(set([round(l, 2) for l in lums_5])) < len(lums_5):
        print("Result: NOT SUITABLE. Duplicate lightness values create indistinguishable colors.\n")
    else:
        print("Result: SUITABLE.\n") # This path won't be taken

    # --- Final Conclusion ---
    print("="*40)
    print("CONCLUSION:")
    if not suitable_plots:
        print("None of the plots are suitable.")
        final_answer = "none"
    else:
        result_str = ",".join(map(str, sorted(suitable_plots)))
        print(f"The only plot interpretable for someone with monochromatic vision is: {result_str}")
        final_answer = result_str
        
    print(f"\n<<<{final_answer}>>>")

if __name__ == "__main__":
    analyze_plots()