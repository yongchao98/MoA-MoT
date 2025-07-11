def analyze_monochromatic_suitability():
    """
    Analyzes color palettes to determine their suitability for monochromatic vision
    by calculating the relative luminance of each color.
    """
    
    print("For a plot to be interpretable by someone with monochromatic vision, the different")
    print("colors used must have distinct lightness values to appear as different shades of gray.")
    print("\nWe will analyze the palette for Plot 2, 'pals::ocean.balance(5)', by calculating")
    print("the luminance of each color. The formula is: Luminance = 0.2126*R + 0.7152*G + 0.0722*B\n")

    # Hex codes for the pals::ocean.balance(5) palette used in Plot 2
    palette_plot_2 = {
        "Color 1": "#296497",
        "Color 2": "#85A6C0",
        "Color 3": "#F6F5F2",
        "Color 4": "#E8B273",
        "Color 5": "#B35A2A"
    }

    luminances = []

    for name, hex_code in palette_plot_2.items():
        # 1. Convert hex to RGB
        hex_val = hex_code.lstrip('#')
        R, G, B = tuple(int(hex_val[i:i+2], 16) for i in (0, 2, 4))

        # 2. Calculate luminance
        luminance = 0.2126 * R + 0.7152 * G + 0.0722 * B
        luminances.append(luminance)

        # 3. Print the equation as requested
        print(f"Analysis for {name} ({hex_code}):")
        print(f"  Luminance = 0.2126 * {R} + 0.7152 * {G} + 0.0722 * {B} = {luminance:.2f}")

    print("\n---\n")
    print("Conclusion for Plot 2:")
    print("The calculated luminance values are approximately: 91.14, 160.86, 245.00, 184.93, 105.46.")
    print("Since these values are all distinct, the colors can be differentiated as shades of gray.")
    print("Plot 2 is suitable.\n")

    print("Analysis of other plots:")
    print("- Plots 1 & 6: Use colors with very similar luminance, making them difficult to distinguish.")
    print("- Plots 3 & 4: Use a constant lightness value (L=60), so all colors appear as the same shade of gray.")
    print("- Plot 5: Uses a set of lightness values that includes a duplicate (20), so two colors will be indistinguishable.")
    
    print("\nTherefore, only Plot 2 uses a color palette suitable for monochromatic vision.")


if __name__ == '__main__':
    analyze_monochromatic_suitability()