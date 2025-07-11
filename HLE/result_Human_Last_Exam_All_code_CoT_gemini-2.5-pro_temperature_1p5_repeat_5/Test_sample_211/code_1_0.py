import math

def solve():
    """
    Calculates and explains the change in the slope of a latitudinal diversity gradient
    after an invasion event, according to Hubbell's unified neutral theory.
    """
    print("This script demonstrates what happens to the slope of a latitudinal diversity gradient")
    print("after a widespread invasive species establishes, based on neutral theory.\n")

    # 1. Define the initial state (pre-invasion) using two representative sites.
    # High-latitude site: Ottawa, Canada
    lat_high = 45.42  # Approximate latitude in degrees North
    div_high_initial = 50   # Hypothetical initial alpha diversity (number of species)

    # Low-latitude site: Tena, Ecuador
    lat_low = -0.99  # Approximate latitude in degrees (South is negative)
    div_low_initial = 250  # Hypothetical initial alpha diversity

    print("--- Initial State (Pre-Invasion) ---")
    print(f"High-Latitude Site (Ottawa): Latitude = {lat_high}°, Initial Diversity = {div_high_initial} species")
    print(f"Low-Latitude Site (Tena): Latitude = {lat_low}°, Initial Diversity = {div_low_initial} species\n")

    # 2. Calculate the initial slope.
    # Slope = (y2 - y1) / (x2 - x1) = (change in diversity) / (change in latitude)
    numerator_initial = div_low_initial - div_high_initial
    denominator_initial = lat_low - lat_high
    slope_initial = numerator_initial / denominator_initial

    print("Calculating the initial slope of the diversity gradient:")
    print(f"Slope = (Diversity at Tena - Diversity at Ottawa) / (Latitude of Tena - Latitude of Ottawa)")
    print(f"Initial Slope Equation: ({div_low_initial} - {div_high_initial}) / ({lat_low} - {lat_high})")
    print(f"Initial Slope = {numerator_initial} / {denominator_initial:.2f} = {slope_initial:.4f}\n")

    # 3. Model the effect of the invasion.
    # Neutral theory predicts a homogenization effect. Diversity is reduced everywhere,
    # leading to a larger absolute species loss in the more diverse tropics.
    # We will model this as a proportional reduction for simplicity.
    reduction_factor = 0.4  # Diversity at all sites is reduced to 40% of its original value.

    div_high_final = int(div_high_initial * reduction_factor)
    div_low_final = int(div_low_initial * reduction_factor)

    print("--- Final State (Post-Invasion) ---")
    print(f"The invasive species reduces diversity at all sites to {int(reduction_factor * 100)}% of the original.")
    print(f"Final High-Latitude Diversity = {div_high_final} species")
    print(f"Final Low-Latitude Diversity = {div_low_final} species\n")

    # 4. Calculate the final slope.
    numerator_final = div_low_final - div_high_final
    denominator_final = lat_low - lat_high # Latitude doesn't change
    slope_final = numerator_final / denominator_final

    print("Calculating the final slope of the diversity gradient:")
    print(f"Final Slope Equation: ({div_low_final} - {div_high_final}) / ({lat_low} - {lat_high})")
    print(f"Final Slope = {numerator_final} / {denominator_final:.2f} = {slope_final:.4f}\n")

    # 5. Conclusion
    print("--- Conclusion ---")
    print(f"The initial slope's magnitude was |{slope_initial:.4f}| = {abs(slope_initial):.4f}.")
    print(f"The final slope's magnitude is |{slope_final:.4f}| = {abs(slope_final):.4f}.")
    print("\nBecause the difference in species diversity between the tropics and temperate zones decreases,")
    print("the slope of the latitudinal diversity gradient becomes shallower (less steep).\n")

solve()
<<<The slope of insect diversity will become shallower (flatter).>>>