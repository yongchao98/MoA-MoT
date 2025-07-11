import hsluv
from colormath.color_objects import sRGBColor
from colormath.color_conversions import convert_color
import pypals

def hex_to_lightness(hex_code):
    """Converts a hex color string to CIE L* (perceptual lightness)."""
    try:
        rgb = sRGBColor.new_from_rgb_hex(hex_code.lstrip('#'))
        lab = convert_color(rgb, 'lab')
        return lab.lab_l
    except Exception:
        return 0

def are_lightness_values_distinguishable(lightness_list, min_diff=5):
    """
    Checks if lightness values are all different by at least a minimum threshold.
    """
    if len(lightness_list) < 2:
        return True
    
    sorted_lightness = sorted(lightness_list)
    for i in range(len(sorted_lightness) - 1):
        if (sorted_lightness[i+1] - sorted_lightness[i]) < min_diff:
            return False
    return True

def analyze_palettes():
    """
    Analyzes the 6 palettes from the R code to determine their suitability
    for monochromatic vision.
    """
    # --- Palette Definitions ---
    # Plot 1 & 6: ggplot2 default
    pal_1_6 = ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']

    # Plot 2: pals::ocean.balance(5)
    rgb_pal2 = pypals.palettes.ocean_balance.colors_5
    pal_2 = [f"#{r:02x}{g:02x}{b:02x}" for r, g, b in rgb_pal2]

    # Plot 3: hsluv with varying hue/saturation, constant lightness=60
    pal_3 = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip([0, 60, 120, 180, 240], [0, 20, 40, 60, 80])]

    # Plot 4: hsluv with varying hue, constant saturation=10, lightness=60
    pal_4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in [0, 60, 120, 180, 240]]

    # Plot 5: hsluv with lightness sampled from (20,50,70,20,90), which includes a duplicate
    pal_5 = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip([0, 72, 144, 216, 288], [20, 50, 70, 90, 20])]
    
    all_palettes = { 1: pal_1_6, 2: pal_2, 3: pal_3, 4: pal_4, 5: pal_5, 6: pal_1_6 }
    
    print("Analysis of Palettes for Monochromatic Vision\n" + "="*45)
    
    interpretable_plots = []
    
    for plot_num, palette_hex in all_palettes.items():
        print(f"\n--- Analyzing Plot {plot_num} ---")
        lightness_values = [hex_to_lightness(color) for color in palette_hex]
        
        # We can reason about the construction of palettes 3, 4, and 5.
        if plot_num in [3, 4]:
            is_good = False
            reason = "Palette was generated with a constant HSLuv lightness, making all colors identical in grayscale."
        elif plot_num == 5:
            is_good = False
            reason = "Palette was generated with sampled lightness values that include a duplicate, making two colors identical in grayscale."
        else:
            is_good = are_lightness_values_distinguishable(lightness_values, min_diff=5)
            if is_good:
                reason = "Lightness values are unique and sufficiently spaced apart."
            else:
                reason = "Some colors have very similar lightness values, making them hard to distinguish."
        
        lightness_str = ", ".join([f"{l:.1f}" for l in lightness_values])
        print(f"Perceptual Lightness (L*) values: [ {lightness_str} ]")
        print(f"Conclusion: {reason}")

        if is_good:
            interpretable_plots.append(str(plot_num))
            
    print("\n" + "="*45)
    
    if not interpretable_plots:
        final_answer = "none"
    else:
        # Get unique plot numbers and sort them
        final_answer = ",".join(sorted(list(set(interpretable_plots))))
        
    print(f"\nFinal Answer: The plot(s) interpretable for someone with monochromatic vision is/are: {final_answer}")
    return final_answer

# Run the analysis and capture the final answer
final_answer_value = analyze_palettes()
print(f"<<<{final_answer_value}>>>")