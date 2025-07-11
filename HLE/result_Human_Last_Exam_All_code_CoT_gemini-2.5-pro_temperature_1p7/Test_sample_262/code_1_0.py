import textwrap

def hex_to_rgb(hex_color):
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def calculate_luminance(hex_color):
    """
    Calculates the relative luminance of a hex color and prints the steps.
    """
    r, g, b = hex_to_rgb(hex_color)
    
    # Normalize RGB values to the range 0-1
    norm_r, norm_g, norm_b = r / 255, g / 255, b / 255
    
    # Apply the formula for relative luminance
    luminance = 0.2126 * norm_r + 0.7152 * norm_g + 0.0722 * norm_b
    
    # Print the detailed calculation
    print(f"Color {hex_color}:")
    print(f"  R,G,B = ({r}, {g}, {b})")
    # Using the final equation with each number as requested
    print(f"  Luminance = 0.2126 * ({r}/255) + 0.7152 * ({g}/255) + 0.0722 * ({b}/255)")
    print(f"  Luminance = 0.2126 * {norm_r:.2f} + 0.7152 * {norm_g:.2f} + 0.0722 * {norm_b:.2f} = {luminance:.3f}")
    
    return luminance

def analyze_palette(name, hex_codes):
    """Analyzes a palette by calculating and printing luminances."""
    print("-" * 40)
    print(f"Analyzing {name}...")
    print("-" * 40)
    
    luminances = [calculate_luminance(color) for color in hex_codes]
    
    print(f"\nSummary for {name}:")
    print(f"Luminance values: {[round(l, 3) for l in luminances]}")
    return luminances

def main():
    intro_text = """
    For a plot to be interpretable for someone with monochromatic vision, the colors used for different categories must have distinct brightness levels (luminance). We can test this by converting the palette colors to their relative luminance, which corresponds to how bright they appear in grayscale. A suitable palette will have a wide and distinct range of luminance values.
    """
    print(textwrap.dedent(intro_text).strip())

    # Palette from Plot 2: pals::ocean.balance(5)
    # This is a diverging palette, which should have good luminance variation.
    pal1 = ['#00429D', '#5E96C3', '#F6F2E2', '#C78D5B', '#93003A']
    lum1 = analyze_palette("Plot 2 Palette (`ocean.balance`)", pal1)
    
    print("\nConclusion: The luminance values are distinct and span a wide range (from 0.139 to 0.948). This palette is INTERPRETABLE.")

    # Palette from Plot 6: scales::hue_pal()(5)
    # This is a qualitative palette with similar luminance, which is bad for grayscale conversion.
    pal5 = ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']
    lum5 = analyze_palette("Plot 6 Palette (`hue_pal`)", pal5)
    
    print("\nConclusion: The luminance values are all very close together (all between 0.562 and 0.598). This palette is NOT INTERPRETABLE.")

    final_explanation = """
    - Plot 1 uses a default palette very similar to Plot 6, making it unsuitable.
    - Plots 3 & 4 use colors with a constant lightness value, making them indistinguishable.
    - Plot 5 uses colors where two have the same lightness value, causing confusion.
    - Only Plot 2 uses a palette with systematically varied lightness.
    
    Therefore, the only suitable plot is Plot 2.
    """
    print(textwrap.dedent(final_explanation))


if __name__ == '__main__':
    main()
