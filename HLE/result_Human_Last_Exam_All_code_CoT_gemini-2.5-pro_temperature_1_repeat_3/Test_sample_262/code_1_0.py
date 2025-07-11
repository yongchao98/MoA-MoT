import numpy as np
import hsluv
import matplotlib.colors as mcolors

def calculate_luminance_std(hex_colors):
    """Calculates the standard deviation of luminance for a list of hex colors."""
    luminances = []
    for hex_color in hex_colors:
        rgb = mcolors.to_rgb(hex_color)
        # Luminance formula (perceptual brightness)
        luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
        luminances.append(luminance)
    return np.std(luminances)

def analyze_palettes():
    """
    Analyzes the palettes from the R code to find which are suitable for
    monochromatic vision (i.e., have high variance in luminance).
    """
    palettes = {}

    # Plot 1 & 6: Default ggplot2 palette (scales::hue_pal)
    # Varies hue, but luminance is not intentionally varied and can be similar for different colors.
    palettes[1] = ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']
    palettes[6] = palettes[1]

    # Plot 2: Diverging palette (pals::ocean.balance)
    # Designed to vary from dark to light to dark, ensuring high luminance variation.
    palettes[2] = ['#00429d', '#73a2c6', '#f1f1f1', '#e67962', '#b2182b']

    # Plot 3: HSLuv with constant lightness (L=60)
    # Perceptually uniform lightness means all colors will have nearly identical luminance.
    hues = [0, 60, 120, 180, 240]
    sats = [h / 3 for h in hues]
    palettes[3] = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip(hues, sats)]
    
    # Plot 4: HSLuv with constant lightness (L=60)
    # Same principle as Plot 3.
    palettes[4] = [hsluv.hsluv_to_hex((h, 10, 60)) for h in hues]

    # Plot 5: HSLuv with varying lightness
    # Lightness is explicitly varied, ensuring different luminances.
    hues = [0, 72, 144, 216, 288]
    lightnesses = [20, 50, 70, 20, 90] # A representative permutation
    palettes[5] = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(hues, lightnesses)]

    interpretable_plots = []
    # A low standard deviation in luminance (<0.1) indicates a poor palette for monochrome vision.
    # Palettes designed with varying luminance will have a significantly higher std dev.
    luminance_threshold = 0.15

    for plot_num, hex_colors in sorted(palettes.items()):
        std_dev = calculate_luminance_std(hex_colors)
        if std_dev > luminance_threshold:
            interpretable_plots.append(str(plot_num))
            
    print(",".join(interpretable_plots))

if __name__ == '__main__':
    analyze_palettes()