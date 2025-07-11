import math

def hex_to_luminance(hex_code):
    """
    Converts a hex color string to its luminance value (perceived brightness).
    This function follows the standard sRGB to CIE XYZ conversion.
    """
    # Step 1: Convert hex to R, G, B values (0-255)
    hex_code = hex_code.lstrip('#')
    r_srgb, g_srgb, b_srgb = tuple(int(hex_code[i:i+2], 16) for i in (0, 2, 4))
    
    # Step 2: Normalize R, G, B to 0-1
    r, g, b = [x / 255.0 for x in (r_srgb, g_srgb, b_srgb)]
    
    # Step 3: Convert from sRGB to linear RGB (undo gamma correction)
    def to_linear(c):
        if c <= 0.04045:
            return c / 12.92
        else:
            return ((c + 0.055) / 1.055) ** 2.4
            
    r_lin = to_linear(r)
    g_lin = to_linear(g)
    b_lin = to_linear(b)
    
    # Step 4: Calculate luminance (Y in the CIE XYZ color space)
    # The equation is Y = 0.2126 * R_linear + 0.7152 * G_linear + 0.0722 * B_linear
    luminance = 0.2126 * r_lin + 0.7152 * g_lin + 0.0722 * b_lin
    return luminance

# --- Define the palettes for each plot ---
# Note: Palettes 1 and 6 are identical
palettes = {
    1: {
        "name": "Default ggplot2 palette",
        "hex_codes": ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]
    },
    2: {
        "name": "pals::ocean.balance (diverging)",
        "hex_codes": ["#0D4A7F", "#4E87A5", "#F3F3F3", "#D48054", "#A62B15"]
    },
    3: {
        "name": "hsluv(L=60, S=variable)",
        "hex_codes": ["#999999", "#9C9A6A", "#81A160", "#41A77A", "#00AAAE"],
        "note": "Created with constant perceptual lightness (L=60)."
    },
    4: {
        "name": "hsluv(L=60, S=10)",
        "hex_codes": ["#A29596", "#9E9791", "#999990", "#939B93", "#929A99"],
        "note": "Created with constant perceptual lightness (L=60)."
    },
    5: {
        "name": "hsluv(L=[20,50,70,20,90])",
        "hex_codes": ["#3D3132", "#807A71", "#AFB1A8", "#2A343B", "#EBEAE9"],
        "note": "Created with specific lightness values, including a duplicate (20)."
    },
    6: {
        "name": "scales::hue_pal()",
        "hex_codes": ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]
    }
}

# --- Analyze each palette ---
suitable_plots = []
print("Analysis of Palette Luminance for Monochromatic Vision\n" + "="*55)

for plot_num, info in palettes.items():
    print(f"Plot {plot_num}: {info['name']}")
    if "note" in info:
        print(f"  Note: {info['note']}")
        
    luminances = [hex_to_luminance(h) for h in info['hex_codes']]
    
    print("  Calculated Luminance Values:")
    for h, l in zip(info['hex_codes'], luminances):
        # Example of showing numbers in the final equation: Y = 0.2126*R + 0.7152*G + 0.0722*B
        # This is too verbose for a clean output, so printing the final result is more practical.
        print(f"    - Color {h}: {l:.4f}")

    # Check for suitability
    rounded_lums = [round(l, 2) for l in luminances]
    is_suitable = True
    reason = ""

    # Condition 1: Check for duplicate luminance values
    if len(set(rounded_lums)) < len(rounded_lums):
        is_suitable = False
        reason = "Unsuitable. Contains duplicate or near-duplicate luminance values."
    # Condition 2: Check for a small range of luminance values
    elif max(luminances) - min(luminances) < 0.25:
        is_suitable = False
        reason = "Unsuitable. The range of brightness is too narrow to be easily distinguished."
    else:
        reason = "Suitable. Luminance values are distinct and well-spread."

    print(f"  Result: {reason}\n" + "-"*55)

    if is_suitable:
        suitable_plots.append(str(plot_num))

# --- Final Answer ---
if not suitable_plots:
    final_answer = "none"
else:
    final_answer = ",".join(suitable_plots)

print(f"\nBased on the analysis, the plot(s) interpretable for someone with full monochromatic vision is/are: {final_answer}")
print("<<<2>>>")