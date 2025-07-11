import numpy as np
from colorspacious import cspace_convert
import hsluv
import pals

def hex_to_srgb1(hex_color):
    """Converts a hex color string to an sRGB tuple with values in [0, 1]."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))

def analyze_palette(plot_num, name, palette_hex, notes=""):
    """Analyzes a palette and prints its suitability for monochromatic vision."""
    pal_srgb1 = [hex_to_srgb1(h) for h in palette_hex]
    pal_lab = cspace_convert(pal_srgb1, "sRGB1", "CIELAB")
    lightness_values = sorted([round(lab[0], 1) for lab in pal_lab])
    
    print(f"--- Plot {plot_num}: {name} ---")
    print(f"Lightness (L*) values: {lightness_values}")

    unique_lightness_count = len(set(lightness_values))
    is_suitable = unique_lightness_count == len(palette_hex)
    
    if not is_suitable:
        print("Result: NOT SUITABLE. Contains colors with identical lightness.")
    elif max(lightness_values) - min(lightness_values) < 40:
        print("Result: NOT SUITABLE. Lightness values are not well-separated.")
    else:
        print("Result: SUITABLE. Lightness values are unique and well-separated.")
        return True
        
    return False

# --- Define Palettes ---
# Plot 1 & 6: scales::hue_pal()(5) -> Uses HCL with constant lightness l=65
pal1_6 = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

# Plot 2: pals::ocean.balance(5) -> A diverging palette
pal2 = pals.ocean.balance_5

# Plot 3: Custom HSLuv with constant lightness l=60
pal3 = [hsluv.hsluv_to_hex((h, h/3, 60)) for h in [0, 60, 120, 180, 240]]

# Plot 4: Custom HSLuv with constant lightness l=60
pal4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in [0, 60, 120, 180, 240]]

# Plot 5: Custom HSLuv with non-unique lightness values [20, 50, 70, 20, 90]
pal5 = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip([0, 72, 144, 216, 288], [20, 50, 70, 20, 90])]

# --- Run Analysis ---
palettes_to_check = [
    (1, "Default ggplot", pal1_6),
    (2, "pals::ocean.balance", pal2),
    (3, "Custom HSLuv (Constant Lightness)", pal3),
    (4, "Custom HSLuv (Constant Lightness)", pal4),
    (5, "Custom HSLuv (Repeated Lightness)", pal5),
    (6, "scales::hue_pal", pal1_6),
]

print("Analyzing palettes for monochromatic vision interpretability...\n")
suitable_plots = []
for plot_num, name, pal in palettes_to_check:
    if analyze_palette(plot_num, name, pal):
        suitable_plots.append(str(plot_num))
    print()

final_answer = ",".join(suitable_plots) if suitable_plots else "none"

print(f"The plot(s) using a color palette interpretable for someone with full monochromatic vision are: {final_answer}")
<<<2>>>