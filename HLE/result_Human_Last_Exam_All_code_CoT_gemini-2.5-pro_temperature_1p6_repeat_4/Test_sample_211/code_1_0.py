import numpy as np

def explain_diversity_change():
    """
    This function models the change in the latitudinal diversity gradient
    following an invasive species event, according to Hubbell's Neutral Theory.
    """

    # Approximate latitudes for the 12 sites from the equator outwards
    # Tena, Cali, BCI, Sarapiqui, Managua, Yoro, Belmopan, Merida, Miami, Charlotte, State College, Ottawa
    latitudes = np.array([0, 3, 9, 10, 12, 15, 17, 21, 26, 35, 41, 45])
    
    # --- Step 1: Model the Initial State ---
    # We assume a strong initial latitudinal gradient.
    # Let's create a hypothetical linear relationship for insect diversity.
    # Diversity = 300 - 5 * Latitude
    initial_diversity = 300 - 5 * latitudes
    
    # Calculate the slope of the initial gradient using linear regression
    initial_slope, initial_intercept = np.polyfit(latitudes, initial_diversity, 1)
    
    print("--- Initial State (Pre-Invasion) ---")
    print(f"A strong latitudinal diversity gradient exists.")
    print(f"Example Diversity at Equator (Tena, ~0 lat): {initial_diversity[0]:.0f} species")
    print(f"Example Diversity at North (Ottawa, ~45 lat): {initial_diversity[-1]:.0f} species")
    print(f"Initial Slope of the Diversity Gradient: {initial_slope:.2f}")
    print("\nInitial Equation: Diversity =", f"{initial_intercept:.2f}", "+", f"({initial_slope:.2f})", "* Latitude\n")


    # --- Step 2: Model the Post-Invasion State ---
    # The invasive species causes biotic homogenization, making communities more similar.
    # This flattens the gradient. We can model this by shrinking the diversity
    # range and making all sites more similar.
    # A simple model: new diversity is the average of the old diversity and the mean diversity.
    mean_diversity = np.mean(initial_diversity)
    final_diversity = (initial_diversity + mean_diversity) / 2
    
    # Calculate the slope of the final gradient
    final_slope, final_intercept = np.polyfit(latitudes, final_diversity, 1)
    
    print("--- Final State (Post-Invasion) ---")
    print("The invasive species causes biotic homogenization, reducing the difference")
    print("in diversity between equatorial and temperate sites.")
    print(f"Example Diversity at Equator (Tena, ~0 lat): {final_diversity[0]:.0f} species")
    print(f"Example Diversity at North (Ottawa, ~45 lat): {final_diversity[-1]:.0f} species")
    print(f"Final Slope of the Diversity Gradient: {final_slope:.2f}")
    print("\nFinal Equation: Diversity =", f"{final_intercept:.2f}", "+", f"({final_slope:.2f})", "* Latitude\n")

    # --- Step 3: Conclusion ---
    print("--- Conclusion ---")
    print("As predicted by neutral theory, biotic homogenization by the invasive species")
    print("causes the magnitude of the slope to decrease, making the gradient flatter (closer to zero).")


explain_diversity_change()

<<<The slope will decrease in magnitude (become flatter, i.e., closer to zero).>>>