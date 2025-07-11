import hsluv
import numpy as np

def get_luminance(hex_code):
    """
    Calculates the luminance of a hex color code.
    The hex code is converted to linear RGB (0-1), then luminance is calculated.
    """
    hex_code = hex_code.lstrip('#')
    r_linear = int(hex_code[0:2], 16) / 255.0
    g_linear = int(hex_code[2:4], 16) / 255.0
    b_linear = int(hex_code[4:6], 16) / 255.0
    
    luminance = 0.2126 * r_linear + 0.7152 * g_linear + 0.0722 * b_linear
    return luminance

def analyze_plots():
    """
    Analyzes the color palettes of the 6 plots to determine their
    suitability for monochromatic vision.
    """
    # Palette definitions corresponding to the plots
    # Note: Plot 1 and 6 use the same default ggplot2 palette
    palettes = {
        1: ("ggplot2 default", ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]),
        2: ("pals::ocean.balance", ["#002D3C", "#0088A0", "#E6E6E6", "#A57321", "#4A2800"]),
        # For plot 3, h=[0,60,120,180,240], s=h/3, l=60
        3: ("hsluv (varying sat)", [hsluv.hsluv_to_hex((h, h/3.0, 60)) for h in np.linspace(0, 240, 5)]),
        # For plot 4, h=[0,60,120,180,240], s=10, l=60
        4: ("hsluv (low sat)", [hsluv.hsluv_to_hex((h, 10, 60)) for h in np.linspace(0, 240, 5)]),
        # For plot 5, h=[0,72,144,216,288], l=sampled from [20,50,70,20,90]
        # The key issue is the duplicate lightness value (20), regardless of sample order.
        5: ("hsluv (random light)", [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(np.arange(0, 289, 72), [20, 50, 70, 20, 90])]),
        6: ("scales::hue_pal()", ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"])
    }

    suitable_plots = []
    print("Analyzing plot palettes for monochromatic interpretability...\n")
    
    for plot_num, (name, hex_codes) in palettes.items():
        luminances = [get_luminance(h) for h in hex_codes]
        
        print(f"--- Plot {plot_num} ({name}) ---")
        for i, (h, l) in enumerate(zip(hex_codes, luminances)):
             print(f"Color {i+1}: {h}, Luminance: {l:.4f}")

        # Criteria for suitability: distinct and spread-out luminance values
        # We check if luminances are nearly identical or have duplicates
        luminances_rounded = [round(l, 3) for l in luminances]
        unique_luminances = set(luminances_rounded)
        luminance_range = max(luminances) - min(luminances)

        if len(unique_luminances) == len(hex_codes) and luminance_range > 0.5:
             print("Result: SUITABLE. Luminance values are distinct and have a wide range.\n")
             suitable_plots.append(str(plot_num))
        elif len(unique_luminances) < len(hex_codes):
             print("Result: UNSUITABLE. Palette has duplicate luminance values.\n")
        elif luminance_range < 0.2:
             print("Result: UNSUITABLE. Luminance values are too close to be distinguished.\n")
        else:
             print("Result: OK, but not ideal. Luminance values are distinct but range is moderate.\n")
             
    if not suitable_plots:
        print("Final Answer: none")
    else:
        print(f"Final Answer: {','.join(suitable_plots)}")

analyze_plots()
<<<2>>>