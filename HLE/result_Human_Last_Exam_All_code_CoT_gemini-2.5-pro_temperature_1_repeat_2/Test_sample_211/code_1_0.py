import numpy as np

def solve():
    """
    Models the effect of an invasive species on the latitudinal diversity gradient
    under the Unified Neutral Theory of Biodiversity.
    """
    # 1. Define the absolute latitudes for the 12 sites, from Tena to Ottawa.
    # The gradient depends on distance from the equator, so we use absolute values.
    latitudes = np.array([
        0.99,  # Tena, Ecuador
        3.45,  # Cali, Colombia
        9.15,  # Barro Colorado Island, Panama
        10.46, # Sarapiqui, Costa Rica
        12.13, # Managua, Nicaragua
        15.14, # Yoro, Honduras
        17.25, # Belmopan, Belize
        20.97, # Merida, Mexico
        25.76, # Miami, Florida
        35.23, # Charlotte, North Carolina
        40.80, # State college, Pennsylvania
        45.42  # Ottawa, Canada
    ])

    # 2. Create a hypothetical initial diversity gradient.
    # We model this as a linear decrease in diversity as latitude increases.
    # Diversity = Intercept - Slope * Latitude
    initial_intercept_param = 160
    initial_slope_param = 2.7
    initial_diversity = initial_intercept_param - initial_slope_param * latitudes

    # 3. Calculate the best-fit slope for the initial diversity gradient.
    # np.polyfit(x, y, 1) returns [slope, intercept] for a linear fit.
    initial_slope, initial_intercept = np.polyfit(latitudes, initial_diversity, 1)

    print("--- Initial State Before Invasion ---")
    print("A diversity gradient exists, with more species near the equator.")
    print("The initial diversity gradient is modeled by the equation:")
    # We print each number in the equation as requested.
    # The slope is negative because diversity decreases as latitude increases.
    print(f"Diversity = {initial_slope:.2f} * Latitude + {initial_intercept:.2f}")

    # 4. Simulate biotic homogenization due to the invasive species.
    # Diversity is reduced at all sites. We model this as diversity dropping
    # to 50% of its original value at each site.
    homogenization_factor = 0.5
    final_diversity = initial_diversity * homogenization_factor

    # 5. Calculate the slope of the final, post-invasion diversity gradient.
    final_slope, final_intercept = np.polyfit(latitudes, final_diversity, 1)

    print("\n--- Final State After Invasion ---")
    print("Biotic homogenization has reduced diversity across all sites.")
    print("The final diversity gradient is modeled by the equation:")
    # We print each number in the final equation.
    print(f"Diversity = {final_slope:.2f} * Latitude + {final_intercept:.2f}")

    print(f"\nThe magnitude of the slope (its steepness) changed from {abs(initial_slope):.2f} to {abs(final_slope):.2f}.")
    print("As shown by the model, the slope has become shallower and is closer to zero.")

solve()