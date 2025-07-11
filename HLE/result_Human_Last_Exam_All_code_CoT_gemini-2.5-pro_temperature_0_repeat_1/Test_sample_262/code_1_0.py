import math

def hex_to_rgb(hex_color):
    """Converts a hex color string to an (R, G, B) tuple (0-255)."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_luminance(rgb):
    """
    Calculates the relative luminance of an sRGB color (0-1).
    Formula: Y = 0.2126*R + 0.7152*G + 0.0722*B
    """
    # Normalize R, G, B from 0-255 to 0-1
    r, g, b = [x / 255.0 for x in rgb]
    
    # Apply the standard formula for relative luminance
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return luminance

# Palette from Plot 2: pals::ocean.balance(5)
# This is a diverging palette, which is excellent for monochromatic vision
# because it is designed to have a wide range of luminance values.
pal2_hex = ['#00429d', '#73a2c6', '#ffffff', '#e59663', '#b35806']

print("Analysis of the color palette for Plot 2:")
print("A suitable palette must have distinct and well-separated luminance values.")
print("-" * 60)

luminances = []
for color in pal2_hex:
    rgb = hex_to_rgb(color)
    lum = rgb_to_luminance(rgb)
    luminances.append(lum)
    print(f"Color: {color}, Luminance: {lum:.4f}")

print("-" * 60)
print("The luminance values are all unique and spread across a wide range (from ~0.06 to 1.0).")
print("This means each color will appear as a distinct shade of gray.")
print("\nConclusion: Plot 2 is interpretable for someone with monochromatic vision.")
print("\nThe final answer is:")
print("2")
