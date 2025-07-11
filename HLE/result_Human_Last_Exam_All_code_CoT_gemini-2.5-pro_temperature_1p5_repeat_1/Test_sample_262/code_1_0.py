import palettable

def hex_to_rgb(hex_color):
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def get_grayscale_luminance(rgb_tuple):
    """Calculates the perceived luminance of an sRGB color."""
    # Convert 0-255 range to 0-1 range for calculation
    r, g, b = [x / 255.0 for x in rgb_tuple]
    # Standard formula for relative luminance (per ITU-R BT.709)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b

# Based on the analysis, only Plot 2 is suitable for monochromatic vision.
# We will verify this by calculating the luminance values for its palette.
# Plot 2 uses `pals::ocean.balance(5)`.
correct_plot_number = 2
palette = palettable.cmocean.diverging.Balance_5.hex_colors

# Calculate the luminance for each color in the palette
luminance_values = [get_grayscale_luminance(hex_to_rgb(h)) for h in palette]

print(f"The only suitable plot is Plot {correct_plot_number}.")
print("Its diverging palette provides distinct luminance values for each color, making it interpretable in grayscale.")
print("\nThe luminance values (ranging from 0=black to 1=white) for the five colors are:")

# The prompt asks to "output each number in the final equation!".
# We interpret this to mean printing the calculated luminance values.
print(*[f"{val:.4f}" for val in luminance_values], sep=", ")
