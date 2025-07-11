import hsluv

def explain_and_solve():
    """
    This function analyzes the color palettes from the R code to determine
    which are suitable for an observer with monochromatic vision (sees in grayscale).

    The plan is as follows:
    1.  Acknowledge that for a plot to be interpretable in monochrome, each color
        representing a different category must have a distinct lightness value.
    2.  Use the HSLuv color space, as its 'L' (Lightness) component is designed
        for perceptual uniformity, making it an excellent proxy for grayscale representation.
    3.  Recreate each of the 6 color palettes described in the R code.
    4.  For each palette, calculate the HSLuv lightness value for every color.
    5.  Check if all lightness values within a single palette are unique.
        - If unique, the palette is suitable.
        - If there are duplicate lightness values, the palette is not suitable because
          different categories would appear as the same shade of gray.
    6.  Report the results for each plot and provide the final answer.
    """
    print("Plan to identify the monochromatically-interpretable plot:")
    print("1. Recreate the 5 distinct color palettes from the R code.")
    print("2. For each palette, determine the perceptual lightness of its colors.")
    print("3. A palette is interpretable if all its lightness values are unique.")
    print("4. Identify the plot(s) that use such a palette.\n")


def check_palette(palette_name, hex_colors):
    """
    Checks if a color palette is suitable for monochromatic vision by
    checking the uniqueness of the HSLuv lightness component.
    Returns True if suitable, False otherwise.
    """
    # The HSLuv Lightness component is a direct measure of perceived lightness.
    lightness_values = [hsluv.hex_to_hsluv(h)[2] for h in hex_colors]
    
    # We round to 2 decimal places to check for effective uniqueness.
    rounded_lightness = [round(l, 2) for l in lightness_values]
    unique_lightness_count = len(set(rounded_lightness))
    
    print(f"--- Analysis for {palette_name} ---")
    print(f"Hex Colors: {hex_colors}")
    print(f"Corresponding HSLuv Lightness Values: {[f'{l:.2f}' for l in lightness_values]}")

    is_interpretable = unique_lightness_count == len(hex_colors)
    if is_interpretable:
        print("Result: INTERPRETABLE in monochrome. All lightness values are distinct.")
    else:
        print("Result: NOT INTERPRETABLE in monochrome. Contains duplicate or near-identical lightness values.")
    print("-" * (len(palette_name) + 21) + "\n")
    return is_interpretable


# --- Main script execution ---

explain_and_solve()

# --- Recreate Palettes ---

# Plot 1 & 6 use the default ggplot palette (`pal5` in the R code)
pal_plot1_and_6 = ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3']

# Plot 2 uses `pal1` from `pals::ocean.balance(5)`
pal_plot2 = ["#3B4992", "#EEF3F7", "#ED7354", "#B4332D", "#6A001A"]

# Plot 3 uses `pal2` where Lightness is a constant 60
hues_p3 = [0, 60, 120, 180, 240]
sats_p3 = [h / 3 for h in hues_p3]
pal_plot3 = [hsluv.hsluv_to_hex((h, s, 60)) for h, s in zip(hues_p3, sats_p3)]

# Plot 4 uses `pal3` where Lightness is a constant 60
hues_p4 = [0, 60, 120, 180, 240]
pal_plot4 = [hsluv.hsluv_to_hex((h, 10, 60)) for h in hues_p4]

# Plot 5 uses `pal4` where lightness values are sampled from a list with a duplicate
# The lightness values `[20, 50, 70, 20, 90]` contain '20' twice, guaranteeing a collision.
hues_p5 = [0, 72, 144, 216, 288]
lightnesses_p5 = [20, 50, 70, 20, 90] 
pal_plot5 = [hsluv.hsluv_to_hex((h, 10, l)) for h, l in zip(hues_p5, lightnesses_p5)]


# --- Analyze each palette and check for suitability ---
results = {}
results[1] = check_palette("Plot 1 (Default ggplot)", pal_plot1_and_6)
results[2] = check_palette("Plot 2 (pals::ocean.balance)", pal_plot2)
results[3] = check_palette("Plot 3 (pal2)", pal_plot3)
results[4] = check_palette("Plot 4 (pal3)", pal_plot4)
results[5] = check_palette("Plot 5 (pal4 with duplicate lightness)", pal_plot5)
results[6] = check_palette("Plot 6 (Default ggplot)", pal_plot1_and_6)

# --- Report the final answer ---
interpretable_plots = [str(k) for k, v in results.items() if v]
final_answer = ",".join(interpretable_plots) if interpretable_plots else "none"

print("\n=================================================")
print("Final Conclusion:")
print(f"The plot number(s) interpretable for someone with monochromatic vision is/are: {final_answer}")
print("=================================================")
