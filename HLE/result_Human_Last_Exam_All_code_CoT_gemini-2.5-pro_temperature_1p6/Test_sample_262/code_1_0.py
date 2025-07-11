import sys

def analyze_palettes():
    """
    Analyzes color palettes to determine their suitability for viewers with
    monochromatic vision.
    A palette is suitable if its colors have distinct luminance values.
    """

    print("Analyzing plot palettes for monochromatic vision suitability.\n")
    print("A palette is suitable if colors have distinct luminance values, making them distinguishable in grayscale.")
    print("Luminance is calculated as: Y = 0.2126*R + 0.7152*G + 0.0722*B\n")

    # Hex codes for palettes where calculation is needed.
    # Palettes are from `scales::hue_pal()(5)` and `pals::ocean.balance(5)` respectively.
    palettes = {
        "Plot 1 and 6": {
            "colors": ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3'],
            "reason": "Default ggplot palette (scales::hue_pal). Designed with constant luminance."
        },
        "Plot 2": {
            "colors": ['#00396A', '#6699CC', '#FFFFFF', '#FFAE73', '#7D3800'],
            "reason": "Diverging palette (pals::ocean.balance), which typically has good luminance variation."
        },
    }

    # Analysis of HSLuv palettes based on their generating code.
    programmatic_palettes = {
        "Plot 3": "Generated with a constant lightness parameter (l=60). The colors will not be distinguishable in grayscale.",
        "Plot 4": "Generated with a constant lightness parameter (l=60). The colors will not be distinguishable in grayscale.",
        "Plot 5": "Generated with lightness values sampled from (20, 50, 70, 20, 90). Two colors will share the same lightness value (l=20), making them indistinguishable."
    }

    suitable_plots = []

    # Analyze palettes by calculating luminance
    for plot_num_str, data in palettes.items():
        print(f"--- Analyzing {plot_num_str} ---")
        print(f"Reasoning: {data['reason']}")
        
        luminances = []
        is_suitable = True

        print("\n  Color   | R  | G  | B  | Luminance Calculation                 | Result")
        print("---------------------------------------------------------------------------------")
        for color in data['colors']:
            r = int(color[1:3], 16)
            g = int(color[3:5], 16)
            b = int(color[5:7], 16)
            
            # Normalize to 0-1 range
            r_lin, g_lin, b_lin = r / 255.0, g / 255.0, b / 255.0

            # Calculate luminance
            luminance = 0.2126 * r_lin + 0.7152 * g_lin + 0.0722 * b_lin
            luminances.append(luminance)
            
            # Output each number in the final equation as requested
            calc_str = f"0.2126*{r_lin:.3f} + 0.7152*{g_lin:.3f} + 0.0722*{b_lin:.3f}"
            print(f"  {color} | {r:3d} | {g:3d} | {b:3d} | {calc_str:<40}| {luminance:.4f}")

        # Check for suitability
        # Check if luminance range is too small
        if max(luminances) - min(luminances) < 0.25:
             is_suitable = False
             print("\nConclusion: NOT SUITABLE. Luminance values are too close to be distinguished.\n")
        else:
             print("\nConclusion: SUITABLE. Luminance values are distinct.\n")
             # "Plot 1 and 6" -> ["1", "6"]
             # "Plot 2" -> ["2"]
             nums = [char for char in plot_num_str if char.isdigit()]
             suitable_plots.extend(nums)
             
    # Print analysis of programmatic palettes
    for plot_num_str, reason in programmatic_palettes.items():
        print(f"--- Analyzing {plot_num_str} ---")
        print(f"Reasoning: {reason}")
        print("Conclusion: NOT SUITABLE.\n")
        
    print("--- FINAL ANSWER ---")
    if not suitable_plots:
        final_answer = "none"
    else:
        final_answer = ",".join(sorted(suitable_plots))
        
    print(f"The plot(s) using a color palette interpretable for someone with full monochromatic vision are: {final_answer}")

if __name__ == '__main__':
    analyze_palettes()