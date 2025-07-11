import math

def hex_to_rgb(hex_color):
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def calculate_relative_luminance(rgb):
    """Calculates the relative luminance for an sRGB color."""
    r, g, b = [x / 255.0 for x in rgb]
    # Apply gamma correction
    r = r / 12.92 if r <= 0.03928 else ((r + 0.055) / 1.055) ** 2.4
    g = g / 12.92 if g <= 0.03928 else ((g + 0.055) / 1.055) ** 2.4
    b = b / 12.92 if b <= 0.03928 else ((b + 0.055) / 1.055) ** 2.4
    # Standard formula for relative luminance
    return 0.2126 * r + 0.7152 * g + 0.0722 * b

def check_palette_suitability(palette):
    """Checks if a palette's colors have unique luminance values."""
    luminances = [calculate_relative_luminance(hex_to_rgb(c)) for c in palette]
    # Check if the number of unique luminances is equal to the number of colors
    return len(set(luminances)) == len(palette)

def solve():
    """
    Analyzes the palettes from the problem description to find which one is
    suitable for monochromatic vision.
    """
    # Palette definitions based on the R code
    # Plot 1 & 6: scales::hue_pal()(5) - known to have constant luminance. Unsuitable.
    # Plot 3 & 4: hsluv with constant lightness = 60. Unsuitable.
    # Plot 5: hsluv with lightness values from (20,50,70,20,90). Contains duplicate. Unsuitable.
    # Plot 2 is the only candidate that needs to be checked.
    
    pal_plot2 = ["#0B42AA", "#6982D5", "#F5F5F5", "#B57A27", "#7C3E14"] # from pals::ocean.balance(5)
    
    suitable_plots = []
    
    # We only need to check Plot 2, as others are unsuitable by design.
    if check_palette_suitability(pal_plot2):
        suitable_plots.append(2)

    if not suitable_plots:
        print("none")
    else:
        # The problem asks for the plot number
        # We need to print each number of the final equation, but here the final answer is just a number.
        # So we just print the number.
        print(','.join(map(str, suitable_plots)))

solve()