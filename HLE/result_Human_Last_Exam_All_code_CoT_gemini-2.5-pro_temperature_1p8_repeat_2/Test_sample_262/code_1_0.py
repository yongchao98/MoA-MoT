def hex_to_luminance(hex_code: str) -> float:
    """
    Converts a hex color string to its perceived luminance.
    Luminance is calculated based on the sRGB and Rec. 709 standards.
    """
    hex_code = hex_code.lstrip('#')
    r, g, b = (int(hex_code[i:i+2], 16) / 255.0 for i in (0, 2, 4))

    def srgb_to_linear(c: float) -> float:
        """Converts an sRGB component to a linear value."""
        if c <= 0.04045:
            return c / 12.92
        else:
            return ((c + 0.055) / 1.055) ** 2.4

    r_lin = srgb_to_linear(r)
    g_lin = srgb_to_linear(g)
    b_lin = srgb_to_linear(b)

    # Perceived luminance calculation
    return 0.2126 * r_lin + 0.7152 * g_lin + 0.0722 * b_lin

def is_palette_monochromatic_friendly(luminances: list[float]) -> bool:
    """
    Checks if a list of luminances is suitable for monochromatic vision.
    It checks for both uniqueness and sufficient range.
    """
    # Rule 1: All luminances must be unique. We round to 2 decimal places
    # to account for very small floating point differences.
    rounded_lums = [round(l, 2) for l in luminances]
    if len(set(rounded_lums)) != len(luminances):
        return False

    # Rule 2: The range of luminances should be large enough to be discernible.
    # A range of at least 0.5 is a good heuristic.
    if (max(luminances) - min(luminances)) < 0.5:
        return False
        
    return True

# Hex codes for each of the 5 distinct palettes used in the 6 plots.
# Palettes generated from the R code logic.
palettes = {
    1: ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3'], # Plot 1: ggplot default
    2: ['#009490', '#56B9AC', '#F0F0F0', '#EDA678', '#DE7043'], # Plot 2: pals::ocean.balance(5)
    3: ['#9e9e9e', '#9b9d75', '#73a563', '#25ab8a', '#00adc8'], # Plot 3: hsluv with varying saturation
    4: ['#9e9d9b', '#9b9e95', '#95a092', '#92a09a', '#959fa0'], # Plot 4: hsluv with constant lightness
    5: ['#3b3b3b', '#85837f', '#bdb9b4', '#393c3e', '#ebebeb'], # Plot 5: hsluv with varying lightness (includes duplicate)
    6: ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']  # Plot 6: Same as Plot 1
}

interpretable_plots = []
for plot_num, palette_hex in palettes.items():
    # Calculate the luminance for each color in the palette
    luminances = [hex_to_luminance(h) for h in palette_hex]
    
    # Check if the set of luminances is friendly to monochromatic vision
    if is_palette_monochromatic_friendly(luminances):
        interpretable_plots.append(str(plot_num))

if not interpretable_plots:
    result = "none"
else:
    result = ",".join(interpretable_plots)

print(result)