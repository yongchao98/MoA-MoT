import hsluv

def hex_to_rgb(hex_str):
    """Converts a hex color string to an (R, G, B) tuple."""
    hex_str = hex_str.lstrip('#')
    return tuple(int(hex_str[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_grayscale_luminance(rgb_tuple):
    """
    Converts an (R, G, B) tuple to a grayscale luminance value (0-255)
    using the standard Rec. 709 formula for perceived brightness.
    """
    r, g, b = rgb_tuple
    # Using the standard formula L = 0.2126*R + 0.7152*G + 0.0722*B
    return 0.2126 * r + 0.7152 * g + 0.0722 * b

def is_palette_monochromatic_friendly(hex_colors, tolerance=5.0):
    """
    Checks if a palette is friendly to monochromatic vision by ensuring
    all colors have distinct grayscale luminance values.
    A tolerance is used to account for very similar but not identical values.
    """
    luminance_values = sorted([rgb_to_grayscale_luminance(hex_to_rgb(c)) for c in hex_colors])
    
    for i in range(len(luminance_values) - 1):
        if abs(luminance_values[i] - luminance_values[i+1]) < tolerance:
            return False
    return True

def main():
    # Recreate the 5 palettes from the R code.
    # The data has 5 groups, so we need 5 colors from each palette.

    # Palette for Plot 1 & 6 (Default ggplot2 palette: scales::hue_pal())
    # This creates colors with the same luminance. We simulate this with HSLuv.
    pal1_6 = [hsluv.hsluv_to_hex((h, 95, 65)) for h in [15, 88, 161, 234, 307]]

    # Palette for Plot 2 (pals::ocean.balance(5))
    # This is a diverging palette, designed for varying lightness.
    pal2 = ['#003c67', '#3c799c', '#e2e2e2', '#b35833', '#7f2800']

    # Palette for Plot 3 (sapply(seq(0, 300, by = 60), \(i) hsluv_hex(i, i/3, 60)))
    # Lightness 'l' is constant at 60.
    hues = [0, 60, 120, 180, 240]
    saturations = [h / 3 for h in hues]
    pal3 = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip(hues, saturations)]
    
    # Palette for Plot 4 (sapply(seq(0, 300, by = 60), \(i) hsluv_hex(i, 10, 60)))
    # Lightness 'l' is constant at 60.
    pal4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in [0, 60, 120, 180, 240]]

    # Palette for Plot 5 (sapply(seq(0, 288, by = 72), \(i) hsluv_hex(i, 10, sample(c(20,50,70,20,90)))))
    # Lightness 'l' values contain a duplicate (20), making two colors indistinguishable in grayscale.
    hues_5 = [0, 72, 144, 216, 288]
    lightness_5 = [20, 50, 70, 20, 90] # The list contains a duplicate lightness value.
    pal5 = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(hues_5, lightness_5)]

    palettes = {
        1: pal1_6,
        2: pal2,
        3: pal3,
        4: pal4,
        5: pal5,
        6: pal1_6
    }
    
    interpretable_plots = []
    for plot_num, pal_hex in palettes.items():
        if is_palette_monochromatic_friendly(pal_hex):
            interpretable_plots.append(str(plot_num))
            
    if interpretable_plots:
        print(",".join(interpretable_plots))
    else:
        print("none")

if __name__ == '__main__':
    main()