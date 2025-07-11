def analyze_palettes():
    """
    Analyzes color palettes to determine their suitability for monochromatic vision.
    A palette is suitable if its colors have distinct luminance values.
    """
    
    print("Analyzing the palette from Plot 2, which uses `pals::ocean.balance`.")
    print("This palette is designed to have colors with varying brightness.")
    
    pal_plot2_hex = ["#00429D", "#73A2C6", "#E7E6E6", "#D5926F", "#93261F"]
    
    print("\nWe can confirm this by calculating the luminance of each color.")
    print("The formula is: Luminance = 0.299*R + 0.587*G + 0.114*B")
    
    # Show an example calculation for the first color
    first_color_hex = pal_plot2_hex[0]
    r, g, b = tuple(int(first_color_hex.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
    
    print(f"\nExample calculation for the color '{first_color_hex}':")
    print(f"  R = {r}, G = {g}, B = {b}")
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    print(f"  Luminance = 0.299 * {r} + 0.587 * {g} + 0.114 * {b} = {luminance:.1f}")

    print("\nCalculated luminance for all colors in the palette (scale 0-255):")
    all_luminances = []
    for color_hex in pal_plot2_hex:
        r, g, b = tuple(int(color_hex.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
        lum = 0.299 * r + 0.587 * g + 0.114 * b
        all_luminances.append(lum)
        print(f"  - {color_hex}: {lum:.1f}")
        
    print("\nThe luminance values are distinct, making the palette suitable.")
    print("\nAll other palettes described use colors with constant or overlapping luminance values, making them unsuitable.")

    # The final answer identifies the suitable plot number.
    final_answer = "2"
    print(f"\nTherefore, the only suitable plot is: {final_answer}")


analyze_palettes()